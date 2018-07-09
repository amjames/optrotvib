import json
from pathlib import Path
import time
from . import cpu_info
from . import psi4_engine
from . import g09_engine

def prepare():
    input_json = json.loads(Path('input.json').read_text())
    output_json = {
            'output': {},
            'raw_output': {},
            'compute_info': {},
            'success': False,
            '_id': input_json.get('_id')
            }
    return input_json, output_json

def wrap_compute_info(f, arg):
    t_start = time.time()
    output = f(arg)
    t_end = time.time()
    output['compute_info'] = {
            'cluster_name': cpu_info.cluster_name(),
            'hostname': cpu_info.hostname(),
            'ncpu': cpu_info.ncore(),
            'memory': "{} MB".format(cpu_info.memory()),
            'walltime': round(t_end - t_start)
            }
    return output


def finish(output_json):
    Path('output.json').write_text(json.dumps(output_json, indent=4))

def main():
    input_json, output_json = prepare()
    r_id = input_json.get('_id')
    if input_json['modelchem']['program'] == 'psi4':
        output_json = wrap_compute_info(psi4_engine.run, input_json)
    if input_json['modelchem']['program'] == 'g09':
        output_json = wrap_compute_info(g09_engine.run, input_json)
    else:
        output_json['raw_output']['error_message'] = "BAD PROG {}".format(input_json['modelchem']['program'])
    output_json['_id'] = r_id
    Path('output.json').write_text(json.dumps(output_json))
