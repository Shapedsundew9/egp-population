"""Test that the population configuration."""
from json import load
from logging import Logger, NullHandler, getLogger
from os import chdir, getcwd, mkdir
from os.path import exists, join
from subprocess import PIPE, CompletedProcess, run
from tempfile import TemporaryDirectory
from typing import Any, Callable, cast
from uuid import uuid4

import pytest
from pypgtable.pypgtable_typing import TableConfigNorm

from egp_population.egp_typing import PopulationConfig, PopulationsConfig
from egp_population.population_config import configure_populations, population_table_default_config, populations_default_config

_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())


# For testing local file population configuration
_FAKE_FITNESS_FUNCTION: str = """
from numpy import single
from typing import Any, Callable


EGP_PROBLEM_CONFIG: dict[str, Any] = {
    "name": "Fake",
    "description": "More fake",
    "inputs": ["str", "str"],
    "outputs": ["float"],
}


def fitness_function(individual: Callable[[int, int], int]) -> single:
    return single(0.0)

"""
_FAKE_SURVIVABILITY_FUNCTION: str = """
from typing import cast
from numpy import bool_, full, single
from numpy.typing import NDArray
from egp_population.population import population
from egp_population.egp_typing import Active, Reference, Survivability, SurvivabilityFunction


def survivability_function(populous: population) -> tuple[Survivability, Active] | tuple[Reference, Survivability, Active]:
    survivability: NDArray[single] = cast(NDArray[single], populous['fitness'])
    active: NDArray[bool_] = full(survivability.shape, False, dtype=bool_)
    return (survivability, active)

"""

# Store the current working directory & change to where Erasmus GP keeps its data
cwd: str = getcwd()

# Create a temporary directory & change to it
temp_dir = TemporaryDirectory()  # pylint: disable=consider-using-with
DIRECTORY_PATH: str = temp_dir.name
chdir(DIRECTORY_PATH)

# Check the verified problem definitions file
PROBLEM_DEFINITIONS: list[dict[str, Any]] = []
problem_definitions_file: str = join(DIRECTORY_PATH, "egp_problems.json")
problem_definitions_file_exists: bool = exists(problem_definitions_file)
if not problem_definitions_file_exists:
    _logger.info(
        f"The egp_problems.json does not exist in {DIRECTORY_PATH}'."
        " Pulling from https://raw.githubusercontent.com/Shapedsundew9/egp-problems/main/egp_problems.json"
    )
    cmdlist: list[str] = ["wget", "-q", "https://raw.githubusercontent.com/Shapedsundew9/egp-problems/main/egp_problems.json"]
    result: CompletedProcess[bytes] = run(cmdlist, stdout=PIPE, check=False, cwd=DIRECTORY_PATH)
    log_level: Callable[..., None] = _logger.error if result.returncode else _logger.info
    log_level(f"\n{result.stdout.decode('utf-8')}")
    if result.returncode:
        raise RuntimeError(f"Command \"{' '.join(cmdlist)}\" FAILED.")

# Load the problems definitions file if it exists
with open(problem_definitions_file, "r", encoding="utf8") as file_ptr:
    PROBLEM_DEFINITIONS = load(file_ptr)
_logger.debug(f"PROBLEM_DEFINITIONS:\n{PROBLEM_DEFINITIONS}")


def test_population_table_config() -> None:
    """Test that the population table configuration is correct."""
    assert isinstance(population_table_default_config(), dict)


def test_configure_populations_empty() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    config["delete_db"] = True
    configure_populations({"configs": []}, PROBLEM_DEFINITIONS, config)


def test_configure_populations_local_file_error() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    config["delete_db"] = True
    pconfig: PopulationConfig = {
        "uid": 1,
        "worker_id": "70c7d47d-d162-48ce-bdd3-7073cd83cde6",
        "created": "2023-09-10T15:36:57.889532Z",
        "egp_problem": DIRECTORY_PATH + "/",
        "name": "Minimal_population_test_configuration",
        "description": "Valid problem path with trailing /.",
        "meta_data": "{'test': 'test'}",
    }
    with pytest.raises(FileExistsError):
        configure_populations({"configs": [pconfig]}, PROBLEM_DEFINITIONS, config)


def test_configure_populations_local_no_name_error() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    config["delete_db"] = True
    pconfig: dict[str, Any] = {
        "uid": 1,
        "worker_id": "70c7d47d-d162-48ce-bdd3-7073cd83cde6",
        "created": "2023-09-10T15:36:57.889532Z",
        "egp_problem": DIRECTORY_PATH + "/",
        "description": "Valid problem path with trailing /.",
        "meta_data": "{'test': 'test'}",
    }
    with pytest.raises(ValueError):
        configure_populations({"configs": [pconfig]}, PROBLEM_DEFINITIONS, config)  # type: ignore


def test_configure_populations_config_exists() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    config["delete_db"] = True
    pconfigs: PopulationsConfig = populations_default_config()
    pconfigs.get("configs", [])[0]["worker_id"] = uuid4()
    _logger.debug(f"pconfigs:\n{pconfigs}")
    p_config_dict, _, __ = configure_populations(pconfigs, PROBLEM_DEFINITIONS, config)

    # Second time through the config will already exist
    pconfigs.get("configs", [])[0]["uid"] = p_config_dict[1]["uid"]
    p_config_dict, _, __ = configure_populations(pconfigs, PROBLEM_DEFINITIONS, config)
    assert isinstance(p_config_dict, dict)
    assert 1 in p_config_dict
    assert p_config_dict[1]["inputs"] == ["int", "int"]


def test_configure_populations_local() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    config["delete_db"] = True
    pconfigs: PopulationsConfig = populations_default_config()
    fake_fitness_function_path: str = "test_module_f/"
    pconfigs.get("configs", [])[0]["worker_id"] = uuid4()
    pconfigs.get("configs", [])[0]["egp_problem"] = fake_fitness_function_path
    mkdir(fake_fitness_function_path)
    with open(join(fake_fitness_function_path, "fitness_function.py"), "w", encoding="utf8") as fileptr:
        fileptr.write(_FAKE_FITNESS_FUNCTION)
    _logger.debug(f"pconfigs:\n{pconfigs}")
    p_config_dict, _, __ = configure_populations(pconfigs, PROBLEM_DEFINITIONS, config)
    assert isinstance(p_config_dict, dict)
    assert 1 in p_config_dict
    assert p_config_dict[1]["inputs"] == ["str", "str"]


def test_configure_populations_local_survivability() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    config["delete_db"] = True
    pconfigs: PopulationsConfig = populations_default_config()
    fake_survivability_function_path: str = "test_module_s/"
    pconfigs.get("configs", [])[0]["worker_id"] = uuid4()
    pconfigs.get("configs", [])[0]["survivability"] = fake_survivability_function_path
    mkdir(fake_survivability_function_path)
    with open(join(fake_survivability_function_path, "survivability_function.py"), "w", encoding="utf8") as fileptr:
        fileptr.write(_FAKE_SURVIVABILITY_FUNCTION)
    _logger.debug(f"pconfigs:\n{pconfigs}")
    p_config_dict, _, __ = configure_populations(pconfigs, PROBLEM_DEFINITIONS, config)
    assert isinstance(p_config_dict, dict)
    assert 1 in p_config_dict
    assert p_config_dict[1]["inputs"] == ["int", "int"]
