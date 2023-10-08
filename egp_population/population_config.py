"""Erasmus GP populations module."""
from __future__ import annotations
import sys
from copy import deepcopy
from json import dumps, load, loads
from logging import DEBUG, Logger, NullHandler, getLogger
from os import chdir, getcwd
from os.path import dirname, join, exists
from subprocess import CompletedProcess, run
from typing import Callable, cast, Any, TYPE_CHECKING

from egp_physics.physics import stablize
from egp_types.eGC import eGC
from egp_types.reference import ref_str
from egp_types.ep_type import interface_definition, ordered_interface_hash, unordered_interface_hash, vtype
from egp_utils.common import default_erasumus_db_config
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


# Circular import & runtime import avoidance
if TYPE_CHECKING:
    from egp_stores.gene_pool import gene_pool


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
    ("inputs", dumps, lambda x: "" if x is None else loads(x)),
    ("outputs", dumps, lambda x: "" if x is None else loads(x)),
)


_POPULATION_TABLE_DEFAULT_CONFIG: TableConfigNorm = raw_table_config_validator.normalized(
    {
        "database": default_erasumus_db_config(),
        "table": "gene_pool_populations",
        "schema": _POPULATION_TABLE_SCHEMA,
        "create_table": True,
        "create_db": True,
        "conversions": _POPULATIONS_CONVERSIONS,
    }
)

_POPULATION_METRICS_TABLE_DEFAULT_CONFIG: TableConfigNorm = raw_table_config_validator.normalized(
    {
        "database": default_erasumus_db_config(),
        "table": "populations_metrics",
        "schema": _POPULATION_METRICS_TABLE_SCHEMA,
        "create_table": True,
        "create_db": True,
    }
)

_POPULATIONS_DEFAULT_CONFIG: PopulationsConfig = {
    "configs": [
        {
            "name": "Number tree problem (default).",
            "description": "This population has been set up using the default configuration.",
            "egp_problem": "b789df5a4654a5165407f84c7509bedca9af8012"
            # 9fec5181e74ccedf96f2a0c664b85b44ef30f331 includes requirements.txt
        }
    ]
}


_FIND_POPULATION_CONFIG_SQL: str = "WHERE {name} = {config_name}"


def population_table_default_config() -> TableConfigNorm:
    """A copy of the populations config pypgtable table configuration."""
    return deepcopy(_POPULATION_TABLE_DEFAULT_CONFIG)


def populations_default_config() -> PopulationsConfig:
    """A copy of the populations config."""
    return deepcopy(_POPULATIONS_DEFAULT_CONFIG)


def run_cmd(cmdlist: list[str]) -> None:
    """Run a command and log the output."""
    result: CompletedProcess[bytes] = run(cmdlist, capture_output=True, check=False)
    log_level: Callable[..., None] = _logger.error if result.returncode else _logger.info
    log_level(f"stdout:\n{result.stdout.decode('utf-8')}")
    log_level(f"Response:\n{result}")
    if result.returncode:
        raise RuntimeError(f"Command \"{' '.join(cmdlist)}\" FAILED.")


def normalize_population_config(p_config: PopulationConfig, problem_file: str, p_table: table) -> PopulationConfigNorm:
    """Normalize a population configuration.

    Args
    ----
    p_config: The population configuration to normalize.
    """
    # Check the problem file exists & format it as a module import
    if not exists(problem_file):
        raise FileExistsError(f"Problem file '{problem_file}' does not exist in {getcwd()}.")
    if not problem_file.endswith(".py"):
        raise ValueError(f"Problem file '{problem_file}' must be a python file.")
    module: str = problem_file[:-3].replace("/", ".").replace("-", "_")

    # Import the problem config (inputs, outputs & creator)
    _logger.info(f"Importing the problem config for population {p_config.get('uid')}.")
    _logger.debug(f"from {module} import EGP_PROBLEM_CONFIG as problem_config_{p_config.get('uid')}")
    sys.path.append(getcwd())
    exec(  # pylint: disable=exec-used
        f"from {module} import EGP_PROBLEM_CONFIG as problem_config_{p_config.get('uid')}", globals(), locals()
    )
    exec(f"p_config['inputs'] = problem_config_{p_config.get('uid')}['inputs']", globals(), locals())  # pylint: disable=exec-used
    exec(f"p_config['outputs'] = problem_config_{p_config.get('uid')}['outputs']", globals(), locals())  # pylint: disable=exec-used
    exec(f"p_config['creator'] = problem_config_{p_config.get('uid')}.get('creator')", globals(), locals())  # pylint: disable=exec-used

    # Calculate interface hashes
    in_eps, _, ins = interface_definition(p_config.get("inputs", []), vtype.EP_TYPE_STR)
    out_eps, _, outs = interface_definition(p_config.get("outputs", []), vtype.EP_TYPE_STR)
    p_config["ordered_interface_hash"] = ordered_interface_hash(in_eps, out_eps, ins, outs)
    p_config["unordered_interface_hash"] = unordered_interface_hash(in_eps, out_eps)

    # Update the population record in the DB with the problem
    literals: dict[str, Any] = {
        "i": p_table.encode_value("inputs", p_config.get("inputs")),
        "o": p_table.encode_value("outputs", p_config.get("outputs")),
        "c": p_config.get("creator"),
        "u": p_config.get("uid"),
    }
    p_table.update("{inputs} = {i}, {outputs} = {o}, {creator} = {c}", "{uid} = {u}", literals)

    # Import the fitness function
    _logger.info(f"Importing the fitness function for population {p_config.get('uid')}.")
    exec(  # pylint: disable=exec-used
        f"from {module} import fitness_function as fitness_function_{p_config.get('uid')}", globals(), locals()
    )
    exec(f"p_config['fitness_function'] = fitness_function_{p_config.get('uid')}", globals(), locals())  # pylint: disable=exec-used

    # Import the survivability function
    _logger.info(f"Importing the survivability function for population {p_config.get('uid')}.")
    if p_config.get("survivability", "elitist") not in SURVIVABILITY_FUNCTIONS:
        survivability_file: str = join(p_config.get("survivability", "./"), "survivability_function.py")
        module = survivability_file[:-3].replace("/", ".")
        import_str: str = f"from {module} import survivability_function as survivability_function_{p_config.get('uid')}\n"
        exec(import_str, globals(), locals())  # pylint: disable=exec-used
        exec(  # pylint: disable=exec-used
            f"p_config['survivability_function'] = survivability_function_{p_config.get('uid')}", globals(), locals()
        )
    else:
        p_config["survivability_function"] = SURVIVABILITY_FUNCTIONS[p_config.get("survivability", "elitist")]
    return cast(PopulationConfigNorm, p_config)


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
    p_table: table = table(raw_table_config_validator.normalized(table_config))
    pm_table_config: TableConfigNorm = deepcopy(_POPULATION_METRICS_TABLE_DEFAULT_CONFIG)
    pm_table_config["table"] = p_table.raw.config["table"] + "_metrics"
    pm_table: table = table(pm_table_config)

    # Final bit is importing fitness & survivability functions
    # Need to make sure the repo exists and is at the right commit first
    problem_dict: dict[str, dict[str, Any]] = {p["git_hash"]: p for p in problem_definitions}
    root_wd: str = getcwd()
    p_configs = []
    for config in populations_config.get("configs", []):
        if config.get("name") is None or not config.get("name"):
            raise ValueError("Population configuration must have a name and it cannot be an empty string. Configuration:\n{config}")

        # Does the configuration already exist?
        p_config: dict[str, Any] = p_table.get(config["name"], {})
        if not p_config:
            _logger.info(f"Population configuration '{config['name']}' does not exist. Creating.")
            p_config_norm: PopulationConfig = population_entry_validator.normalized(config)
            assert p_config_norm is not None, population_entry_validator.error_str()
            p_table.insert([p_config_norm], ("uid",), "tuple")
            p_config.update(p_config_norm)
            p_config["uid"] = p_table[config["name"]]["uid"]
        else:
            _logger.info(f"Population configuration '{config['name']}' exists.")

        # Get the problem definition
        _logger.info(f"Loading population name: {config.get('name')}, UID: {config.get('uid')}")
        pdef: dict[str, Any] = problem_dict.get(p_config.get("egp_problem", "./"), {})
        if not pdef:
            _logger.info(f"Problem definition '{p_config.get('egp_problem', './')}' not found in problem definitions.")
            problem_file = join(p_config.get("egp_problem", "./"), "fitness_function.py")
        else:
            # It is a verified git repo. Clone it to 'git_folder' if it is not already.
            # NOTE: git_folder is the name of the repo with '-' replaced with '_' to allow it to be treated as a python module path
            pdef["git_folder"] = pdef["git_repo"].replace("-", "_")
            if not exists(pdef["git_folder"]):
                _logger.info(f"Git repo '{pdef['git_repo']}' does not exist. Cloning.")
                url: str = pdef["git_url"] + ("/", "")[pdef["git_url"].endswith("/")] + pdef["git_repo"] + ".git"
                run_cmd(["git", "clone", "-n", "--depth=1", "--filter=tree:0", url, pdef["git_folder"]])

            # Checkout the problem at the correct commit if needed
            if not exists(join(pdef["git_folder"], pdef["root_path"])):
                _logger.info(f"Git repo '{pdef['git_repo']}' checking out {pdef['root_path']} at commit {pdef['git_hash']}.")
                chdir(pdef["git_folder"])
                run_cmd(["git", "sparse-checkout", ("set", "add")[exists(".git/info/sparse-checkout")], "--no-cone", pdef["root_path"]])
                run_cmd(["git", "checkout", pdef["git_hash"], "--", pdef["root_path"]])
                if not exists(join(pdef["root_path"], "fitness_function.py")):
                    raise FileNotFoundError(f"Problem file '{join(pdef['root_path'], 'fitness_function.py')}' does not exist.")
                chdir(root_wd)

            # Install the requirements if needed
            chdir(join(pdef["git_folder"], pdef["root_path"]))
            if exists("requirements.txt"):
                _logger.info("Installing requirements.")
                run_cmd(["pip", "install", "-r", "requirements.txt"])
            else:
                _logger.info("No requirements.txt file found.")
            chdir(root_wd)
            problem_file: str = join(pdef["git_folder"], pdef["root_path"], "fitness_function.py")

        p_configs.append(normalize_population_config(cast(PopulationConfig, p_config), problem_file, p_table))

    _logger.info(f"{len(p_configs)} Population configurations normalized.")
    return {c["uid"]: cast(PopulationConfigNorm, c) for c in p_configs}, p_table, pm_table


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
    # TODO: Need to add modes for creating new populations
    #     Mixture of subsets of from a set of populations with matching interfaces of the following:
    #         a. Fittest individuals
    #         b. Random individuals
    #         c. Fittest-diverse individuals
    #         d. Specific individuals (any population)
    #         e. eGC's with same interfaces
    #         f. Most survivable individuals
    #     pGC's with the same criteria as above
    num_to_create: int = num if num is not None else population_config["size"]
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
        gpool.pool[rgc["ref"]]["active"] = True
        gpool.pool.update(fgc_dict)
        _logger.debug(f'Created GGCs to add to Gene Pool: {[ref_str(ggc["ref"]) for ggc in (rgc, *fgc_dict.values())]}')
