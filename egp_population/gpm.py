"""The Gene Pool Module."""

from logging import getLogger, NullHandler, DEBUG
from egp_physics.ep_type import asstr


_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)


def name_func(ref):
    ref_str = ref.to_bytes(8,'big', signed=True).hex()
    return f'ref_{ref_str}'

def write_arg(iab, c):
    return "(" + ", ".join([str(c[arg[1]]) if arg[0] == 'C' else arg[0].lower() + "[" + str(arg[1]) + "]" for arg in iab]) + ",)"

def callable_string(gc):
    # TODO: Add a depth parameter that generates functions that have less codons as single functions.
    string = "# ref: " + str(gc['ref']) + "\n"
    string += "# i = (" + ", ".join((asstr(i) for i in gc['igraph'].input_if())) + ")\n"
    string += "def " + name_func(gc['ref'])
    string += "(i):\n" if len(gc['inputs']) else "():\n"
    if _LOG_DEBUG and len(gc['inputs']): string += f"\t_logger.debug(f'{gc['ref']}: i = {{i}}')\n"
    graph = gc['igraph'].app_graph
    if not 'meta_data' in gc or not 'function' in gc['meta_data']:
        c = graph['C'] if 'C' in graph else tuple()
        if gc.get('gca_ref', None) is not None:
            string += "\ta = " + name_func(gc['gca_ref'])
            string += "(" + write_arg(graph['A'], c) + ")\n" if 'A' in gc['graph'] else "()\n"
            if _LOG_DEBUG: string += f"\t_logger.debug(f'{gc['ref']}: a = {{a}}')\n"
        if gc.get('gcb_ref', None) is not None:
            string += "\tb = " + name_func(gc['gcb_ref'])
            string += "(" + write_arg(graph['B'], c) + ")\n" if 'B' in gc['graph'] else "()\n"
            if _LOG_DEBUG: string += f"\t_logger.debug(f'{gc['ref']}: b = {{b}}')\n"
        retval = write_arg(graph['O'], c)
        if _LOG_DEBUG: string += f"\t_logger.debug(f'{gc['ref']}: return {retval} = {{{retval}}}')\n"
        string += "\treturn " + retval + "\n\n\n"
    else:
        format_dict = {'c' + str(i): v for i, v in enumerate(graph['C'])} if 'C' in graph else {}
        format_dict.update({'i' + str(i): 'i[{}]'.format(i) for i in range(len(gc['inputs']))})
        code = gc['meta_data']['function']['python3']['0']
        if 'code' in code: string += "\t" + code['code'].format(**format_dict) + "\n"
        formated_inline = code['inline'].format(**format_dict)
        if _LOG_DEBUG: string += f"\t_logger.debug(f'{gc['ref']}: {formated_inline} = {{{formated_inline}}}')\n"
        string += "\treturn (" + formated_inline + ",)\n\n\n"
    if _LOG_DEBUG: _logger.debug(f"Callable string created:\n{string}")
    return string

def create_callable(gc):
    # Import imports into the global namespace
    if 'meta_data' in gc and 'function' in gc['meta_data']:
        python = gc['meta_data']['function']['python3']['0']
        if 'imports' in python:
            for impt in python['imports']:
                if impt['name'] not in globals():
                    string = "from {module} import {object} as {name}\n".format(**impt)
                    if _LOG_DEBUG: _logger.debug(f"New imports executable: {string}")
                    exec(string, globals())

    exec(callable_string(gc), globals())
    gc['exec'] = globals()[name_func(gc['ref'])]
