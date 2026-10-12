"""Shared fixtures. Tests that need a table skip, with the command that builds it, if it is absent."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glyco import schema  # noqa: E402
from glyco.cohorts import HUMAN, MOUSE, PDAC_TMA  # noqa: E402

SOURCE = ROOT / "glyco"


def _require(path):
    if not path.exists():
        pytest.skip(f"{path.name} not found; run `python -m glyco extract` (needs raw_images/)")


@pytest.fixture(scope="session")
def wide():
    _require(HUMAN.table("wide"))
    return schema.read_wide(HUMAN)


@pytest.fixture(scope="session")
def long():
    _require(HUMAN.table("long"))
    return schema.read_long(HUMAN)


@pytest.fixture(scope="session")
def cases():
    _require(HUMAN.table("cases"))
    return schema.read_cases(HUMAN)


@pytest.fixture(scope="session")
def mouse_wide():
    _require(MOUSE.table("wide"))
    return schema.read_wide(MOUSE)


@pytest.fixture(scope="session")
def mouse_cases():
    _require(MOUSE.table("cases"))
    return schema.read_cases(MOUSE)


@pytest.fixture(scope="session")
def tma():
    _require(PDAC_TMA.table("wide"))
    return schema.read_wide(PDAC_TMA)


@pytest.fixture(scope="session")
def source_files():
    return sorted(SOURCE.rglob("*.py"))
