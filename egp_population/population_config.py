"""Erasmus GP populations module."""
from copy import deepcopy
from json import dumps, load, loads
from logging import DEBUG, Logger, NullHandler, getLogger
from os import chdir, getcwd, listdir
from os.path import dirname, join, exists
from subprocess import PIPE, CompletedProcess, run
from typing import LiteralString, Callable, cast, Any

from egp_physics.physics import stablize
from egp_stores.gene_pool import gene_pool
from egp_types.eGC import eGC
from egp_types.ep_type import vtype
from egp_types.reference import ref_str
from pypgtable import table
from pypgtable.pypgtable_typing import (
    Conversions,
    TableConfig,
    TableConfigNorm,
    TableSchema,
)
from pypgtable.validators import raw_table_config_validator

from .population_validator import population_entry_validator
from .egp_typing import PopulationConfig, PopulationsConfig, PopulationConfigNorm
from .survivability import SURVIVABILITY_FUNCTIONS


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
# TODO: Add a _LOG_CONSISTENCY which additionally does consistency checking
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


with open(
    join(dirname(__file__), "formats/population_table_format.json"),
    "r",
    encoding="utf8",
) as file_ptr:
    _POPULATION_TABLE_SCHEMA: TableSchema = load(file_ptr)
with open(
    join(dirname(__file__), "formats/population_metrics_format.json"),
    "r",
    encoding="utf8",
) as file_ptr:
    _POPULATION_METRICS_TABLE_SCHEMA: TableSchema = load(file_ptr)


_POPULATIONS_CONVERSIONS: Conversions = (
    ("inputs", dumps, loads),
    ("outputs", dumps, loads),
)


_POPULATION_TABLE_DEFAULT_CONFIG: TableConfigNorm = raw_table_config_validator.normalized(
    {
        "database": {"dbname": "erasmus"},
        "table": "gene_pool_populations",
        "schema": _POPULATION_TABLE_SCHEMA,
        "create_table": True,
        "create_db": True,
        "conversions": _POPULATIONS_CONVERSIONS,
    }
)

_POPULATION_METRICS_TABLE_DEFAULT_CONFIG: TableConfigNorm = raw_table_config_validator.normalized(
    {
        "database": {"dbname": "erasmus"},
        "table": "populations_metrics",
        "schema": _POPULATION_METRICS_TABLE_SCHEMA,
        "create_table": True,
        "create_db": True,
    }
)

_POPULATIONS_DEFAULT_CONFIG: PopulationsConfig = {
    "configs": [
        {
            "inputs": [],
            "outputs": [],
        }
    ]
}


_POPULATION_CONFIGS_SQL: LiteralString = "WHERE {population_hash} = ANY({hashes})"
_POPULATION_CONFIG_EXISTS_SQL: LiteralString = "WHERE {population_hash} = ANY({hashes})"
_POPULATION_UID_EXISTS_SQL: LiteralString = "WHERE {uid} = ANY({uids})"


def population_table_default_config() -> TableConfigNorm:
    """A copy of the populations config pypgtable table configuration."""
    return deepcopy(_POPULATION_TABLE_DEFAULT_CONFIG)


def populations_default_config() -> PopulationsConfig:
    """A copy of the populations config."""
    return deepcopy(_POPULATIONS_DEFAULT_CONFIG)


def run_cmd(cmdlist: list[str]) -> None:
    """Run a command and log the output."""
    result: CompletedProcess[bytes] = run(cmdlist, stdout=PIPE, check=False)
    log_level: Callable[..., None] = _logger.error if not result.returncode else _logger.info
    log_level(f"\n{result.stdout.decode('utf-8')}")
    if result.returncode:
        raise RuntimeError(f"Command \"{' '.join(cmdlist)}\" FAILED.")


def configure_populations(
    populations_config: PopulationsConfig,
    problem_definitions: list[dict[str, Any]],
    table_config: TableConfig | TableConfigNorm = _POPULATION_TABLE_DEFAULT_CONFIG,
) -> tuple[dict[int, PopulationConfigNorm], table, table]:
    """Configure populations for runtime.

    Args
    ----
    populations_config: List of minimal population configurations & general parameters.
    table_config: The population pypgtable table configuration.
    """
    p_table_config: TableConfigNorm = raw_table_config_validator.normalized(table_config)
    p_table: table = table(p_table_config)
    pm_table_config: TableConfigNorm = deepcopy(_POPULATION_METRICS_TABLE_DEFAULT_CONFIG)
    pm_table_config["table"] = p_table_config["table"] + "_metrics"
    pm_table: table = table(pm_table_config)

    # If a population config has a UID defined then it should already exist
    # Find it and update the config.
    pcs: list[PopulationConfig] = populations_config.get("configs", [])
    defined_uids: dict[int, PopulationConfig] = {config["uid"]: config for config in pcs if "uid" in config}
    for defined_config in p_table.select(_POPULATION_UID_EXISTS_SQL, {"uids": list(defined_uids.keys())}):
        defined_uids[defined_config["uid"]].update(defined_config)

    # Its possible that a UID was set in populations_config that does not exist.
    # If enough config was defined to find or create a new population then the UID will be updated.
    # If there is not enough the validator below will error.
    for config in populations_config.get('configs', []):
        if not population_entry_validator.validate(config):
            assert False, f"Population configuration is not valid:\n{population_entry_validator.error_str()}"

    # To be fully normalized the configs must get existing UIDs from the DB or have the DB create them
    # Population configurations must have the worker_id populated. If the population config
    # has already been created by another worker then the worker_id will be over written with the correct UUID.
    configs: tuple[PopulationConfigNorm, ...] = tuple(population_entry_validator.normalized(p) for p in populations_config.get('configs', []))
    hashes: list[bytes] = [p["population_hash"] for p in configs if "population_hash" in p]
    existing: set[bytes] = set()
    if hashes:
        existing = set(
            p_table.select(
                _POPULATION_CONFIG_EXISTS_SQL,
                {"hashes": hashes},
                ("population_hash",),
                "tuple",
            )
        )
    new_configs = tuple(config for config in configs if config.get("population_hash", 0) in set(hashes) - existing)
    if new_configs:
        p_table.insert(new_configs)
        p_configs: dict[int, PopulationConfig] = {c["uid"]: c for c in p_table.select(_POPULATION_CONFIGS_SQL, {"hashes": hashes})}
    else:
        p_configs = {}

    # Final bit is importing fitness & survivability functions
    # Need to make sure the repo exists and is at the right commit first
    problem_dict: dict[str, dict[str, Any]] = {p['git_hash']: p for p in problem_definitions}
    root_wd: str = getcwd()
    for config in p_configs.values():
        _logger.info(f"Loading population name: {config.get('Name')}, UID: {config.get('uid')}, hash: {config.get('population_hash')}")

        # Get the problem definition
        pdef: dict[str, Any] = problem_dict.get(config.get("egp_problem", "./"), {})
        if not pdef:
            problem_file = join(config.get("egp_problem", "./"), "fitness_function.py")
        else:
            # It is a verified git repo. Clone it if it is not already.
            if not exists(pdef["git_repo"]):
                _logger.info(f"Git repo '{pdef['git_repo']}' does not exist. Cloning.")
                url: str = pdef["git_url"] + ("/", "")[pdef["git_url"].endswith("/")] + pdef["git_repo"] + ".git"
                run_cmd(["git", "clone", "-n", "--depth=1", "--filter=tree:0", url])

            # Checkout the problem at the correct commit if needed
            if not exists(join(pdef["git_repo"], pdef["root_path"])):
                _logger.info(f"Git repo '{pdef['git_repo']}' checking out {pdef['root_path']} at commit {pdef['git_hash']}.")
                chdir(pdef["git_repo"])
                operation: str = "add" if listdir(".") else "set"
                run_cmd(["git", "sparse-checkout", operation, "--no-cone", pdef['root_path']])
                run_cmd(["git", "checkout", pdef['git_hash'], "--", pdef['root_path']])
                if not exists(join(pdef["root_path"], "fitness_function.py")):
                    raise FileNotFoundError(f"Problem file '{join(pdef['root_path'], 'fitness_function.py')}' does not exist.")
                chdir(root_wd)

            # Install the requirements if needed
            chdir(join(pdef["git_repo"], pdef["root_path"]))
            _logger.info("Installing requirements.")
            if exists("requirements.txt"):
                run_cmd(["pip", "install", "-r", "requirements.txt"])
            chdir(root_wd)

            problem_file: str = join(pdef["git_repo"], pdef["root_path"], "fitness_function.py")

        # Check the problem file exists & format it as a module import
        if not exists(problem_file):
            raise FileExistsError(f"Problem file '{problem_file}' does not exist.")
        if not problem_file.endswith(".py"):
            raise ValueError(f"Problem file '{problem_file}' must be a python file.")
        module: str = problem_file[:-3].replace("/", ".")


        # Import the problem config
        _logger.info(f"Importing the problem config for population {config.get('uid')}.")
        exec(f"from {module} import EGP_PROBLEM_CONFIG as problem_config_{config.get('uid')}", globals())  # pylint: disable=exec-used
        exec(f"config.update(problem_config_{config.get('uid')})")  # pylint: disable=exec-used

        # Import the fitness function
        _logger.info(f"Importing the fitness function for population {config.get('uid')}.")
        exec(f"from {module} import fitness_function as fitness_function_{config.get('uid')}", globals())  # pylint: disable=exec-used
        exec(f"config['fitness_function'] = fitness_function_{config.get('uid')}")  # pylint: disable=exec-used

        # Import the survivability function
        _logger.info(f"Importing the survivability function for population {config.get('uid')}.")
        if config.get("survivability", "elitist") not in SURVIVABILITY_FUNCTIONS:
            survivability_file: str = join(config.get("survivability", "./"), "survivability_function.py")
            module = survivability_file[:-3].replace("/", ".")
            import_str: str = f"from {module} import survivability_function as survivability_function_{config.get('uid')}\n"
            exec(import_str, globals())  # pylint: disable=exec-used
            exec(f"config['survivability_function'] = survivability_function_{config.get('uid')}")  # pylint: disable=exec-used
        else:
            config["survivability_function"] = SURVIVABILITY_FUNCTIONS[config.get("survivability", "elitist")]

    _logger.info(f"{len(p_configs)} Population configurations normalized.")
    return {u: cast(PopulationConfigNorm, c) for u, c in p_configs.items()}, p_table, pm_table


def new_population(population_config: PopulationConfigNorm, gpool: gene_pool, num: int | None = None) -> None:
    """Create num target population GC's.

    Args
    ----
    population_config: Population definition.
    glib: The geneomic library to get GC's for stablizing the new population individuals.
    gpool: The gene pool to store the new individuals into.
    num: The number of GC's to create. If None it is population_config['active_size']
    """
    # Check the population index exists
    num_to_create: int = num if num is not None else population_config.get("active_size", 100)
    _logger.info(f"Creating {num_to_create} GC's for population index: {population_config['uid']}")

    # Create some new gGC's & mark them as individuals too.
    # This may require pulling new agc's from the genomic library through steady state exceptions
    # in the stabilise() in which case we need to pull in all dependents not already in the
    # gene pool.
    for _ in range(num_to_create):
        egc: eGC = eGC(
            inputs=population_config["inputs"],
            outputs=population_config["outputs"],
            vt=vtype.EP_TYPE_STR,
        )
        rgc, fgc_dict = stablize(gpool, egc)

        rgc["population_uid"] = population_config["uid"]
        gpool.pool[rgc["ref"]] = rgc
        gpool.pool.update(fgc_dict)
        _logger.debug(f'Created GGCs to add to Gene Pool: {[ref_str(ggc["ref"]) for ggc in (rgc, *fgc_dict.values())]}')
