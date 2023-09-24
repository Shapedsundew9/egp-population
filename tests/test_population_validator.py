"""Unit tests for the population_validator module."""
from logging import DEBUG, Logger, NullHandler, getLogger
from json import load
from typing import Any
from os.path import dirname, join
import pytest

from egp_population.population_validator import population_entry_validator


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


with open(join(dirname(__file__), "data/population_config_data.json"), "r", encoding="utf-8") as fileptr:
    _TEST_DATA: list[dict[str, Any]] = load(fileptr)


def test_init() -> None:
    """Test that the module is imported correctly."""
    assert population_entry_validator is not None


def test_population_validator_validation() -> None:
    """Test that the population validator validates the minimal population entry correctly."""
    if not population_entry_validator.validate(_TEST_DATA[0]):
        assert False, f"{population_entry_validator.error_str()}"
    assert True


def test_population_validator_normalization_0() -> None:
    """Test that the population validator normalizes the minimal population entry correctly."""
    config = population_entry_validator.normalized(_TEST_DATA[0])
    if not population_entry_validator.validate(config):
        assert False, f"{population_entry_validator.error_str()}"
    assert True


def test_population_validator_normalization_1() -> None:
    """Test that the population validator normalizes the minimal population entry correctly."""
    config = population_entry_validator.normalized(_TEST_DATA[1])
    if not population_entry_validator.validate(config):
        assert False, f"{population_entry_validator.error_str()}"
    assert True


@pytest.mark.parametrize("i", range(2, 3))
def test_population_validator_valid_cases(i: int) -> None:
    """Test that the population validator correctly validates."""
    if not population_entry_validator.validate(_TEST_DATA[i]):
        assert False, f"{population_entry_validator.error_str()}"
    assert True


@pytest.mark.parametrize("i", range(3, len(_TEST_DATA)))
def test_population_validator_error_cases(i: int) -> None:
    """Test that the population validator detects error in the confguration."""
    assert not population_entry_validator.validate(_TEST_DATA[i])
