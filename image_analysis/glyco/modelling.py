"""Cross-validation, the fitted classifiers, and the origin prediction.

Every cross-validation leaves out one specimen at a time. The unit that is held out is the
`specimen_id` column, and what that means differs between cohorts for a structural reason:

* Human cohort: a patient contributes two or three fields of view of a single section, all
  with the same tissue, staining run, imaging session and label. Holding out one view would
  leave its siblings in training, so a view-grouped estimate would measure recognition of a
  section the model had already seen. The patient is the unit.
* Mouse cohort: each animal contributed one section per organ, so no two sections of a class
  come from the same animal. The section is the unit.

One model object per analysis is serialised, and it is the only source of projections and
predictions: no later step fits its own.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import joblib
import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    cohen_kappa_score,
    confusion_matrix,
)
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict

from .cohorts import NEGATIVE_CONTROL
from .config import GLYCANS, MODEL_DIR
from .selection import features, load

#: Serialised classifiers, by analysis.
MODEL_FILES = {
    "primary": "lda_primary.joblib",
    "subtype": "lda_bc_subtype.joblib",
    "phenotype": "lda_pdac_phenotype.joblib",
    "mouse": "lda_mouse_organ.joblib",
}

TITLES = {
    "primary": "Nine-class cancer classification",
    "subtype": "Breast cancer subtype",
    "phenotype": "PDAC glycophenotype",
    "mouse": "Murine organ classification",
}

#: The analyses on the human cohort (the mouse cohort is handled in `mouse.py`).
HUMAN_ANALYSES = ("primary", "subtype", "phenotype")


def new_estimator() -> LinearDiscriminantAnalysis:
    """The classifier, defined in one place.

    Singular-value solver, class priors taken from the data, no feature scaling: the five
    channels already share an intensity scale, and scaling would move the discriminant
    coordinates the projections are drawn from without changing any classification.
    """
    return LinearDiscriminantAnalysis(solver="svd")


# --------------------------------------------------------------------------------------
# Cross-validation
# --------------------------------------------------------------------------------------


@dataclass
class CVResult:
    """Out-of-fold predictions for one analysis, plus everything derived from them."""

    analysis: str
    labels: list[str]
    y_true: np.ndarray
    y_pred: np.ndarray
    groups: np.ndarray
    n_specimens: int
    n_views: int
    warnings: list[str] = field(default_factory=list)

    @property
    def accuracy(self) -> float:
        return float(accuracy_score(self.y_true, self.y_pred))

    @property
    def kappa(self) -> float:
        return float(cohen_kappa_score(self.y_true, self.y_pred))

    @property
    def n_folds(self) -> int:
        return len(set(self.groups))

    def confusion(self, *, normalise: bool = False) -> pd.DataFrame:
        """Pooled confusion matrix; rows are the reference, columns the prediction."""
        cm = confusion_matrix(self.y_true, self.y_pred, labels=self.labels)
        if normalise:
            cm = cm / cm.sum(axis=1, keepdims=True) * 100
        return pd.DataFrame(cm, index=self.labels, columns=self.labels).rename_axis(
            index="reference", columns="prediction"
        )

    def per_class(self) -> pd.DataFrame:
        """Precision, recall (= sensitivity), F1, support and specificity for each class.

        Specificity is one-versus-rest: of the patches that are not of the class, the share
        that were not assigned to it.
        """
        report = classification_report(
            self.y_true, self.y_pred, labels=self.labels, output_dict=True, zero_division=0
        )
        table = pd.DataFrame({k: report[k] for k in self.labels}).T.rename_axis("label")
        cm = confusion_matrix(self.y_true, self.y_pred, labels=self.labels)
        true_positive = np.diag(cm)
        false_positive = cm.sum(axis=0) - true_positive
        true_negative = cm.sum() - cm.sum(axis=0) - cm.sum(axis=1) + true_positive
        table["specificity"] = true_negative / (true_negative + false_positive)
        return table

    def per_specimen(self) -> pd.Series:
        """Out-of-fold accuracy for each held-out specimen."""
        return (
            pd.DataFrame({"y": self.y_true, "p": self.y_pred, "g": self.groups})
            .groupby("g")
            .apply(lambda d: (d.y == d.p).mean(), include_groups=False)
            .rename("accuracy")
            .rename_axis("specimen_id")
        )

    def per_case(self) -> pd.DataFrame:
        """One row per held-out specimen: its class, patch counts and out-of-fold accuracy."""
        frame = pd.DataFrame({"specimen_id": self.groups, "label": self.y_true,
                              "correct": self.y_true == self.y_pred})
        table = frame.groupby("specimen_id", sort=True).agg(
            label=("label", "first"), n_patches=("correct", "size"), n_correct=("correct", "sum"))
        table["accuracy"] = self.per_specimen()
        return table.reset_index()

    def spread(self) -> dict[str, float]:
        """The distribution of per-specimen accuracy, which is the account of variability."""
        acc = self.per_specimen()
        return {
            "median": float(acc.median()),
            "q1": float(acc.quantile(0.25)),
            "q3": float(acc.quantile(0.75)),
            "min": float(acc.min()),
            "max": float(acc.max()),
            "n_below_half": int((acc < 0.5).sum()),
        }

    def summary(self) -> dict:
        """One row for the cross-validation summary table."""
        return {
            "analysis": self.analysis,
            "n_specimens": self.n_specimens,
            "n_views": self.n_views,
            "n_patches": len(self.y_true),
            "n_classes": len(self.labels),
            "n_folds": self.n_folds,
            "accuracy": self.accuracy,
            "kappa": self.kappa,
            **{f"specimen_accuracy_{k}": v for k, v in self.spread().items()},
        }


def cross_validate_frame(df: pd.DataFrame, analysis: str) -> CVResult:
    """Leave-one-specimen-out cross-validation of the patches in `df`.

    Flags any class represented by fewer than three specimens: with so few, one held-out
    specimen moves that class's recall by a large step, and the point estimate should not be
    read on its own. The out-of-fold predictions are pooled.
    """
    X, y, groups = features(df)

    warnings = []
    per_class = df.groupby("label").specimen_id.nunique()
    thin = per_class[per_class < 3]
    if len(thin):
        warnings.append(
            f"{', '.join(f'{k} has only {v} patients' for k, v in thin.items())} — "
            f"one held-out patient moves that class's recall by up to "
            f"{', '.join(f'{100 / v:.0f} points' for v in thin)}."
        )

    pred = cross_val_predict(new_estimator(), X, y, groups=groups, cv=LeaveOneGroupOut())
    return CVResult(
        analysis=analysis,
        labels=sorted(pd.unique(y)),
        y_true=y,
        y_pred=pred,
        groups=groups,
        n_specimens=df.specimen_id.nunique(),
        n_views=df.view_id.nunique(),
        warnings=warnings,
    )


def cross_validate(analysis: str, *, patches=None) -> CVResult:
    """Leave-one-patient-out cross-validation of one human analysis."""
    return cross_validate_frame(load(analysis, patches=patches), analysis)


# --------------------------------------------------------------------------------------
# Final fit
# --------------------------------------------------------------------------------------


def fit_bundle(df: pd.DataFrame, analysis: str) -> dict:
    """Fit on every specimen in `df` and return the model with what a consumer needs.

    The bundle carries the estimator, the feature order, the class order and the specimens
    it saw, so nothing downstream has to re-derive them.
    """
    X, y, _ = features(df)
    estimator = new_estimator().fit(X, y)
    return {
        "analysis": analysis,
        "estimator": estimator,
        "features": list(GLYCANS),
        "classes": list(estimator.classes_),
        "training_specimens": sorted(df.specimen_id.unique()),
        "n_patches": len(df),
        "explained_variance_ratio": list(estimator.explained_variance_ratio_),
    }


def save_model(bundle: dict) -> None:
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, MODEL_DIR / MODEL_FILES[bundle["analysis"]])


def load_model(analysis: str = "primary") -> dict:
    """Read back a serialised classifier."""
    return joblib.load(MODEL_DIR / MODEL_FILES[analysis])


def fit_final(analysis: str, *, patches=None, save: bool = True) -> dict:
    """Fit one human analysis on all of its specimens and serialise it."""
    bundle = fit_bundle(load(analysis, patches=patches), analysis)
    if save:
        save_model(bundle)
    return bundle


def project(bundle: dict, df: pd.DataFrame) -> pd.DataFrame:
    """Discriminant coordinates LD1..LDn for the given patches."""
    scores = bundle["estimator"].transform(df[bundle["features"]].to_numpy(dtype="float64"))
    return pd.DataFrame(
        scores, columns=[f"LD{i + 1}" for i in range(scores.shape[1])], index=df.index
    )


# --------------------------------------------------------------------------------------
# Metastatic origin
# --------------------------------------------------------------------------------------


def predict_origin(bundle: dict, *, patches=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score the metastatic lesions against the nine-class classifier.

    These specimens are external by construction: the training selection asks for tumour
    tissue, so no fold and no fit ever saw them. Returns per-patch predictions and a summary
    in which every class appears, so that a zero is visible rather than absent.
    """
    df = load("origin", patches=patches)
    assert not set(df.specimen_id) & set(bundle["training_specimens"]), "metastasis reached training"

    predicted = bundle["estimator"].predict(df[bundle["features"]].to_numpy(dtype="float64"))
    per_patch = df[["specimen_id", "label", "view_id"]].assign(predicted=predicted)

    counts = per_patch.groupby(["label", "predicted"]).size().rename("patches").reset_index()
    counts["percent"] = counts.patches / counts.groupby("label").patches.transform("sum") * 100
    everything = pd.MultiIndex.from_product(
        [sorted(per_patch.label.unique()), list(bundle["classes"])], names=["label", "predicted"]
    )
    summary = (
        counts.set_index(["label", "predicted"])
        .reindex(everything, fill_value=0)
        .reset_index()
        .sort_values(["label", "patches"], ascending=[True, False])
    )
    return per_patch, summary


def negative_control(summary: pd.DataFrame) -> pd.DataFrame:
    """Patches assigned to the organ the lesions were resected from, which should be none."""
    return (
        summary[summary.predicted == NEGATIVE_CONTROL]
        .rename(columns={"label": "origin"})[["origin", "predicted", "patches", "percent"]]
        .reset_index(drop=True)
    )


# --------------------------------------------------------------------------------------


def run_human(patches: pd.DataFrame, *, save: bool = True) -> dict:
    """Every human-cohort result in one dict: cross-validations, final models, origin."""
    results: dict = {"cv": {}, "final": {}}
    for analysis in HUMAN_ANALYSES:
        results["cv"][analysis] = cross_validate(analysis, patches=patches)
        results["final"][analysis] = fit_final(analysis, patches=patches, save=save)

    per_patch, summary = predict_origin(results["final"]["primary"], patches=patches)
    results["origin"] = {
        "per_patch": per_patch,
        "summary": summary,
        "negative_control": negative_control(summary),
    }
    return results
