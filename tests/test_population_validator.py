"""Unit tests for the population_validator module."""
from json import load
from typing import Any
from os.path import dirname, join
import pytest

from egp_population.population_validator import population_entry_validator


with open(join(dirname(__file__), "data/population_config_data.json"), "r", encoding="utf-8") as fileptr:
    _TEST_DATA: list[dict[str, Any]] = load(fileptr)


def test_init() -> None:
    """Test that the module is imported correctly."""
    assert population_entry_validator is not None


def test_population_validator_validation() -> None:
    """Test that the population validator validates the minimal population entry correctly."""
    assert population_entry_validator.validate(_TEST_DATA[0])


def test_population_validator_normalization_0() -> None:
    """Test that the population validator normalizes the minimal population entry correctly."""
    config = population_entry_validator.normalized(_TEST_DATA[0])
    del config['ordered_interface_hash']
    del config['population_hash']
    assert population_entry_validator.validate(config)


def test_population_validator_normalization_1() -> None:
    """Test that the population validator normalizes the minimal population entry correctly."""
    config = population_entry_validator.normalized(_TEST_DATA[1])
    del config['ordered_interface_hash']
    del config['population_hash']
    assert population_entry_validator.validate(config)


@pytest.mark.parametrize("i", range(2, len(_TEST_DATA)))
def test_population_validator_error_cases(i: int) -> None:
    """Test that the population validator detects error in the confguration."""
    assert not population_entry_validator.validate(_TEST_DATA[i])
