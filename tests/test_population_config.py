"""Test that the population configuration."""
from egp_population.population_config import new_population, configure_populations, population_table_config


def test_population_table_config() -> None:
    """Test that the population table configuration is correct."""
    assert isinstance(population_table_config(), dict)

