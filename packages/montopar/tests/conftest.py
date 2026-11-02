"""Shared test fixtures for montopar."""
import json
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"


def load_fixture(name: str) -> dict:
    """Load a JSON fixture file."""
    with open(FIXTURES / name, encoding="utf-8") as f:
        return json.load(f)
