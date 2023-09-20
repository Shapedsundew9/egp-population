"""Test that the population configuration."""
from egp_population.population_config import new_population, configure_populations, population_table_default_config
from pypgtable.pypgtable_typing import TableConfigNorm 


def test_population_table_config() -> None:
    """Test that the population table configuration is correct."""
    assert isinstance(population_table_default_config(), dict)


def test_configure_populations() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    configure_populations(config)