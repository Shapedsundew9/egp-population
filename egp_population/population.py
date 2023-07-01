"""Erasmus GP populations module."""
from copy import deepcopy
from json import dumps, load, loads
from logging import DEBUG, Logger, NullHandler, getLogger
from os import chdir, getcwd, makedirs
from os.path import dirname, isdir, join, isfile
from subprocess import PIPE, CompletedProcess, run
from typing import LiteralString, Callable

from egp_physics.physics import stablize
from egp_stores.gene_pool import gene_pool
from egp_stores.genomic_library import genomic_library
from egp_types.eGC import eGC
from egp_types.ep_type import vtype
from egp_types.reference import ref_str
from pypgtable import table
from pypgtable.typing import Conversions, TableConfig, TableConfigNorm, TableSchema
from pypgtable.validators import raw_table_config_validator

from .population_validator import population_entry_validator
from .egp_typing import PopulationConfig, PopulationsConfig, PopulationConfigNorm, cast_pc_to_pcn


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
# TODO: Add a _LOG_CONSISTENCY which additionally does consistency checking
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


with open(join(dirname(__file__), "formats/population_table_format.json"), "r", encoding="utf8") as file_ptr:
    _POPULATION_TABLE_SCHEMA: TableSchema = load(file_ptr)
with open(join(dirname(__file__), "formats/population_metrics_table_format.json"), "r", encoding="utf8") as file_ptr:
    _POPULATION_METRICS_TABLE_SCHEMA: TableSchema = load(file_ptr)


_POPULATIONS_CONVERSIONS: Conversions = (
    ('inputs', dumps, loads),
    ('outputs', dumps, loads)
)


_DEFAULT_POPULATIONS_CONFIG: TableConfigNorm = raw_table_config_validator.normalized({
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool_populations',
    'schema': _POPULATION_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
    'conversions': _POPULATIONS_CONVERSIONS
})

_DEFAULT_POPULATION_METRICS_CONFIG: TableConfigNorm = raw_table_config_validator.normalized({
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'populations_metrics',
    'schema': _POPULATION_METRICS_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
})


_POPULATION_CONFIGS_SQL: LiteralString = '{population_hash} = ANY({hashes})'
_POPULATION_CONFIG_EXITS_SQL: LiteralString = '{population_hash} = ANY({config_hashes})'
_POPULATION_UID_EXISTS_SQL: LiteralString = '{uid} = ANY({uids})'


def population_table_config() -> TableConfigNorm:
    """A copy of the populations config pypgtable table configuration."""
    return deepcopy(_DEFAULT_POPULATIONS_CONFIG)


def configure_populations(populations_config: PopulationsConfig,
                          table_config: TableConfig | TableConfigNorm = _DEFAULT_POPULATIONS_CONFIG
                          ) -> tuple[dict[int, PopulationConfigNorm], table, table]:
    """Configure populations for runtime.

    Args
    ----
    populations_config: List of minimal population configurations & general parameters.
    table_config: The population pypgtable table configuration.
    """
    p_table_config: TableConfigNorm = raw_table_config_validator.normalized(table_config)
    p_table: table = table(p_table_config)
    pm_table_config: TableConfigNorm = deepcopy(_DEFAULT_POPULATION_METRICS_CONFIG)
    pm_table_config['table'] = p_table_config['table'] + '_metrics'
    pm_table: table = table(pm_table_config)

    # If a population config has a UID defined then it should already exist
    # Find it and update the config.
    pcs: list[PopulationConfig] = populations_config.get('configs', [])
    defined_uids: dict[int, PopulationConfig] = {config['uid']: config for config in pcs if 'uid' in config}
    for defined_config in p_table.select(_POPULATION_UID_EXISTS_SQL, {'uids': list(defined_uids.keys())}):
        defined_uids[defined_config['uid']].update(defined_config)

    # Its possible that a UID was set in populations_config that does not exist.
    # If enough config was defined to find or create a new population then the UID will be updated.
    # If there is not enough the validator below will error.
    for config in populations_config:
        if not population_entry_validator.validate(config):
            assert False, f"Population configuration is not valid:\n{population_entry_validator.error_str()}"

    # To be fully normalized the configs must get existing UIDs from the DB or have the DB create them
    # Population configurations must have the worker_id populated. If the population config
    # has already been created by another worker then the worker_id will be over written with the correct UUID.
    configs: tuple[PopulationConfig, ...] = tuple(population_entry_validator.normalized(p) for p in populations_config)
    hashes: list[bytes] = [p['population_hash'] for p in configs]
    existing: set[bytes] = set(p_table.select(_POPULATION_CONFIG_EXITS_SQL, {'hashes': hashes}, ('population_hash',), 'tuple'))
    p_table.insert((config for config in configs if config['population_hash'] in set(hashes) - existing))
    p_configs: dict[int, PopulationConfig] = {c['uid']: c for c in p_table.select(_POPULATION_CONFIGS_SQL, {'hashes': hashes})}

    # Final bit is importing fitness & survivability functions
    # Need to make sure the repo exists and is at the right commit first
    for config in p_configs.values():
        _logger.info(f"Loading population name: {config.get('Name')}, UID: {config.get('uid')}, hash: {config.get('population_hash')}")
        if config.get('git_repo') is not None:  # Then all git fields are not None.
            cwd: str = getcwd()
            git_repo: str = config.get('git_repo', '.')
            git_hash: str = config.get('git_hash', '')
            git_path: str = config.get('git_path', '.')

            # Get to the right path
            if not isdir(git_path):
                _logger.info(f"'{git_path}' does not exist. Creating.")
                makedirs(git_path)
            _logger.info(f"Population definition in git repo '{git_repo}' at commit {git_hash}")
            _logger.info(f"Looking for git repo '{git_repo}' in '{git_path}'")
            chdir(git_path)

            # If the repo exists check it is what we are expecting
            if isdir(config.get('git_repo', '__invalid_dir__')):
                _logger.info(f"Git repo folder '{git_repo}' exists. Using current checkout.")
                chdir(config.get('git_repo', '.'))
                result: CompletedProcess[bytes] = run(['git', 'rev-parse', '--verify', 'HEAD'], stdout=PIPE, check=False)
                log_level: Callable[..., None] = _logger.error if not result.returncode else _logger.info
                log_level('\n' + result.stdout.decode('utf-8'))  # pylint: disable=W1201
                assert result.returncode == 0, (f"Is {git_repo} a valid git repo? ",
                                                "'git rev-parse --verify HEAD' produced a non-zero return code.")
                hash_str: str = result.stdout.decode('utf-8')
                assert len(hash_str) == 41 or len(hash_str) == 65, "Unexpected length of returned hash '{hash_str}'."
                hash_str = hash_str[:-1]
                if populations_config.get('error_on_commit_hash_mismatch', True):
                    assert hash_str == git_hash, f"Git has is '{hash_str}' was expecting \'{git_hash}\'."
                    _logger.info("Git checkout hash matches expectations.")
                elif hash_str != git_hash:
                    _logger.warning(f"Git has is '{hash_str}' was expecting \'{git_hash}\'.")
                else:
                    _logger.info("Git checkout hash matches expectations.")
                _logger.info("Not installing requirements of population git repo as repo already exists.")
            else:
                # Clone the repo
                url: str = config.get('git_url', '') + '/' + git_repo + '.git'
                result: CompletedProcess[bytes] = run(['git', 'clone', url], stdout=PIPE, check=False)
                log_level: Callable[..., None] = _logger.error if not result.returncode else _logger.info
                log_level('\n' + result.stdout.decode('utf-8'))  # pylint: disable=W1201
                assert result.returncode == 0, (f"Is {url} a valid git repo URL? ",
                                                f"'git clone {url}' produced a non-zero return code.")

                # Checkout the commit
                chdir(git_repo)
                result = run(['git', 'checkout', git_hash], stdout=PIPE, check=False)
                log_level: Callable[..., None] = _logger.error if not result.returncode else _logger.info
                log_level('\n' + result.stdout.decode('utf-8'))  # pylint: disable=W1201
                assert result.returncode == 0, (f"Is {git_hash} a valid commit hash for this repo? ",
                                                f"'git checkout {git_hash}' produced a non-zero return code.")

                # Install the dependencies
                if isfile('requirements.txt'):
                    result: CompletedProcess[bytes] = run(['pip', 'install', '-r', 'requirements.txt'], stdout=PIPE, check=False)
                    log_level: Callable[..., None] = _logger.error if not result.returncode else _logger.info
                    log_level('\n' + result.stdout.decode('utf-8'))  # pylint: disable=W1201
                    assert result.returncode == 0, ("Non-zero return code trying to install python requirements. ",
                                                    "'pip install -r requirements.txt'")
                else:
                    _logger.info(f"No requirements.txt file exists for population {config.get('uid')}.")
            chdir(cwd)

        # Import the preload function
        _logger.info(f"Importing the preload function for population {config.get('uid')}.")
        preload_import: str = config.get('preload_import', '__configuration_error__')
        import_str: str = f"from {preload_import} import preload_function as preload_function_{config.get('uid')}\n"
        exec(import_str, globals())
        _logger.info(f"Executing preload function for population {config.get('uid')}.")
        exec(f"preload_function_{config.get('uid')}()")

        # Import the fitness function
        _logger.info(f"Importing the fitness function for population {config.get('uid')}.")
        fitness_import: str = config.get('fitness_import', '__configuration_error__')
        import_str: str = f"from {fitness_import} import fitness_function as fitness_function_{config.get('uid')}\n"
        exec(import_str, globals())
        exec(f"config['fitness_function'] = fitness_function_{config.get('uid')}")

        # Import the survivability function
        _logger.info(f"Importing the survivability function for population {config.get('uid')}.")
        survivability_import: str = config.get('survivability_import', '__configuration_error__')
        import_str = f"from {survivability_import} import survivability_function as survivability_function_{config.get('uid')}\n"
        exec(import_str, globals())
        exec(f"config['survivability_function'] = survivability_function_{config.get('uid')}")

    _logger.info(f'{len(p_configs)} Population configurations normalized.')
    return {u: cast_pc_to_pcn(c) for u, c in p_configs.items()}, p_table, pm_table


def new_population(population_config: PopulationConfigNorm, glib: genomic_library, gpool: gene_pool, num: int | None = None) -> None:
    """Create num target population GC's.

    Args
    ----
    population_config: Population definition.
    glib: The geneomic library to get GC's for stablizing the new population individuals.
    gpool: The gene pool to store the new individuals into.
    num: The number of GC's to create. If None it is population_config['size']
    """
    # Check the population index exists
    num_to_create: int = num if num is not None else population_config['size']
    _logger.info(f"Creating {num_to_create} GC's for population index: {population_config['uid']}")

    # Create some new gGC's & mark them as individuals too.
    # This may require pulling new agc's from the genomic library through steady state exceptions
    # in the stabilise() in which case we need to pull in all dependents not already in the
    # gene pool.
    for _ in range(num_to_create):
        retry_count: int = 0
        rgc: dict | None = None
        fgc_dict: dict = {}
        while rgc is None:
            if retry_count:
                _logger.info(f'eGC random creation failed. Retrying...{retry_count}')
            retry_count += 1
            egc: eGC = eGC(inputs=population_config['inputs'], outputs=population_config['outputs'], vt=vtype.EP_TYPE_STR)
            rgc, fgc_dict = stablize(glib, egc)
            if retry_count == 3:
                raise ValueError(f"Failed to create eGC with inputs = {population_config['inputs']} and outputs"
                                 f" = {population_config['outputs']} {retry_count} times in a row.")

        rgc['population_uid'] = population_config['uid']
        gpool.pool[rgc['ref']] = rgc
        gpool.pool.update(fgc_dict)
        _logger.debug(f'Created GGCs to add to Gene Pool: {[ref_str(ggc["ref"]) for ggc in (rgc, *fgc_dict.values())]}')
