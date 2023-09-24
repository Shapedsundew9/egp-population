"""Test that the population configuration."""
from json import load
from logging import Logger, NullHandler, getLogger
from os import chdir, getcwd
from os.path import exists, join
from tempfile import TemporaryDirectory
from typing import Any, Callable
from subprocess import run, PIPE, CompletedProcess

from pypgtable.pypgtable_typing import TableConfigNorm
from requests import Response, get

from egp_population.egp_typing import PopulationConfig
from egp_population.population_config import (configure_populations,
                                              population_table_default_config)

_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())


# Store the current working directory & change to where Erasmus GP keeps its data
cwd: str = getcwd()

# Create a temporary directory & change to it
temp_dir = TemporaryDirectory()
DIRECTORY_PATH: str = temp_dir.name
chdir(DIRECTORY_PATH)

# Check the verified problem definitions file
PROBLEM_DEFINITIONS: list[dict[str, Any]] = []
problem_definitions_file: str = join(DIRECTORY_PATH, 'egp_problems.json')
problem_definitions_file_exists: bool = exists(problem_definitions_file)
if not problem_definitions_file_exists:
    _logger.info(f"The egp_problems.json does not exist in {DIRECTORY_PATH}'."
                 " Pulling from https://github.com/Shapedsundew9/egp-problems/blob/main/egp_problems.json")
    cmdlist: list[str] = ["wget", "-q", "https://raw.githubusercontent.com/Shapedsundew9/egp-problems/main/egp_problems.json"]
    result: CompletedProcess[bytes] = run(cmdlist, stdout=PIPE, check=False)
    log_level: Callable[..., None] = _logger.error if not result.returncode else _logger.info
    log_level(f"\n{result.stdout.decode('utf-8')}")
    if result.returncode:
        raise RuntimeError(f"Command \"{' '.join(cmdlist)}\" FAILED.")

# Load the problems definitions file if it exists
if problem_definitions_file_exists:
    with open(problem_definitions_file, "r", encoding="utf8") as file_ptr:
        PROBLEM_DEFINITIONS = load(file_ptr)


def test_population_table_config() -> None:
    """Test that the population table configuration is correct."""
    assert isinstance(population_table_default_config(), dict)


def test_configure_populations_empty() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    configure_populations({"configs": []}, PROBLEM_DEFINITIONS, config)


def test_configure_populations() -> None:
    """Test that the population configuration is correct."""
    config: TableConfigNorm = population_table_default_config()
    pconfig: PopulationConfig = {
        "uid": 1,
        "worker_id": "70c7d47d-d162-48ce-bdd3-7073cd83cde6",
        "created": "2023-09-10T15:36:57.889532Z",
        "egp_problem": DIRECTORY_PATH + "/",
        "name": "Minimal_population_test_configuration",
        "description": "Valid problem path with trailing /.",
        "meta_data": "{'test': 'test'}"
    }
    configure_populations({"configs": [pconfig]}, PROBLEM_DEFINITIONS, config)
