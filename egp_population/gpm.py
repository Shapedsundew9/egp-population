"""The Gene Pool Module."""

def name_func(ref):
    ref_str = ref.to_bytes(8,'big', signed=True).hex()
    return f'ref_{ref_str}'

def write_arg(iab, c):
    return "(" + ", ".join([str(c[arg[1]]) if arg[0] == 'C' else arg[0].lower() + "[" + str(arg[1]) + "]" for arg in iab]) + ",)"

def callable_string(gc):
    # TODO: Add a depth parameter that generates functions that have less codons as single functions.
    string = "# ref: " + str(gc['ref']) + "\ndef " + name_func(gc['ref']) + "(i):\n"
    graph = gc['igraph'].app_graph
    if not 'meta_data' in gc or not 'function' in gc['meta_data']:
        c = graph['C'] if 'C' in graph else tuple()
        if gc.get('gca_ref', None) is not None and 'A' in graph:
            string += "\ta = " + name_func(gc['gca_ref']) + "(" + write_arg(graph['A'], c) + ")\n"
        if gc.get('gcb_ref', None) is not None and 'B' in graph:
            string += "\tb = " + name_func(gc['gcb_ref']) + "(" + write_arg(graph['B'], c) + ")\n"
        string += "\treturn " + write_arg(graph['O'], c) + "\n\n\n"
    else:
        format_dict = {'c' + str(i): v for i, v in enumerate(graph['C'])} if 'C' in graph else {}
        format_dict.update({'i' + str(i): 'i[{}]'.format(i) for i in range(len(gc['inputs']))})
        code = gc['meta_data']['function']['python3']['0']
        if 'code' in code: string += "\t" + code['code'].format(**format_dict) + "\n"
        string += "\treturn (" + code['inline'].format(**format_dict) + ",)\n\n\n"
    return string

def create_callable(gc):
    # Import imports into the global namespace
    if 'meta_data' in gc and 'function' in gc['meta_data']:
        python = gc['meta_data']['function']['python3']['0']
        if 'imports' in python:
            for impt in python['imports']:
                impt['name'] = impt['module'].replace('.', '_') + '_' + impt['object'].replace('.', '_')
                if impt['name'] not in globals(): exec("from {module} import {object} as {name}\n".format(**impt), globals())

    exec(callable_string(gc), globals())
    gc['exec'] = globals()[name_func(gc['ref'])]
