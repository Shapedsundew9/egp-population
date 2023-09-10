"""Direct imports."""
from .population import population
from .population_validator import population_entry_validator
from .population_config import population_table_config, new_population, configure_populations

__all__: list[str] = ["population", "population_entry_validator", "population_table_config", "new_population", "configure_populations"]
