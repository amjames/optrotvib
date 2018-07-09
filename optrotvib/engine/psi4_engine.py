import sys
import os
import traceback

from pathlib import Path

from . import cpu_info
from victor.api import Molecule as vicMol
from victor import util

def extract_rotations(all_vars_dict):
    rots_keys = [k for k in all_vars_dict.keys() if "SPECIFIC ROTATION" in k]
    rots = []
    for k in rots_keys:
        rot = {}
        value = all_vars_dict[k]
        method, specific, rotation, p_gauge, at, wl_unit = k.split()
        gauge = p_gauge.strip(')').strip('(')
        wl = int(wl_unit.strip('NM'))
        rot['value'] = value
        rot['wavelength'] = wl
        rot['gauge'] = gauge
        rots.append(rot)
    return rots


def run(input_json):
    # detect custom psi4, if we can
    psiapi_path = os.environ.get('PSIAPI_PATH')
    if psiapi_path:
        sys.path.insert(1, psiapi_path)

    import psi4

    mol_json = input_json.pop('molecule')
    mc_json = input_json.pop('modelchem')
    output_json = {}
    output_json['raw_output'] = {}
    output_json['success'] = False
    output_json['output'] = {}

    # set some psi4 global stuff
    # Scratch files
    scr_path = os.environ.get("WORK")
    psi4_io = psi4.core.IOManager.shared_object()
    psi4_io.set_default_path(scr_path)

    # raw output file
    outfile_path = Path('output.dat')
    psi4.core.set_output_file(str(outfile_path), False)

    # node (or configured) memory and cores (after setting output file so the changes are registered there?)
    psi4.set_memory(str(cpu_info.memory()) + " MB")
    psi4.set_num_threads(cpu_info.ncore())

    # set the basis
    basis = mc_json.get('basis')
    psi4.core.set_global_option('basis', basis)

    # set any other options
    for k, v in mc_json.get('program_options', {}).items():
        psi4.core.set_global_option(k, v)

    try:
        # create the molecule
        psi_mol = psi4.geometry(vicMol(mol_json, dtype='json').to_string())
        # extract name
        method = mc_json.get('method')
        driver = mc_json.get('driver')
        if driver == 'gradient':
            ret, wfn = psi4.gradient(method, return_wfn=True, molecule =psi_mol)
        elif driver == 'hessian':
            ret, wfn = psi4.hessian(method, return_wfn = True, molecule = psi_mol)
        elif driver == 'rotation':
            ret, wfn = psi4.properties(method, return_wfn = True, molecule = psi_mol, properties=['rotation'])
        else:
            raise RuntimeError("invalid driver {}: Validation should have caught this".format(driver))

        grad = wfn.gradient()
        hess = wfn.hessian()
        rotations = extract_rotations(psi4.core.get_variables())
        if hess:
            output_json['output']['hessian'] = util.pack_json_field(hess.to_array())
        if grad:
            output_json['output']['gradient'] = util.pack_json_field(grad.to_array())
        if rotations:
            output_json['output']['rotations'] = rotations
        output_json['output']['all_variables'] = psi4.core.get_variables()
        output_json['success'] = True
        output_json['raw_output']['outfile'] = outfile_path.read_text()
        return output_json
    except Exception as e:
        output_json['success'] = False
        output_json['raw_output']['error'] = "\n".join(traceback.format_exception(*sys.exc_info()))
        return output_json
