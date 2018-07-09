import sys
import re
import os
from pathlib import Path
import subprocess as sp
import numpy as np
import traceback

from . import cpu_info


from victor.constants import physconst
from victor import util

# Collection helpers
def _array_from_fchk(fchk_text, array_name):
    matcher = re.compile(r'\A(?P<title>{})\s+R\s+N=\s+(?P<nele>\d+)\Z'.format(array_name))
    fchk_lines = fchk_text.split('\n')
    start_line = 0
    nline = 0
    found_match = False
    for i, line in enumerate(fchk_lines):
        match = matcher.match(line)
        if match is not None:
            found_match = True
            start_line = i +1
            nline = int(match.group('nele'))//5 + (1 * bool(int(match.group('nele'))%5))
    if found_match:
        data = np.array([float(x) for x in " ".join(fchk_lines[start_line:start_line+nline]).split()])
        return data
    else:
        return None

def compute_rotation(w_au, rot_tensor, mw):
    hbar = physconst['h'] / (2*np.pi)
    prefactor = -72e6 * (hbar**2) * physconst['na'] / physconst['c']**2 / physconst['me']**2
    return prefactor * (w_au**2) * np.trace(rot_tensor) / mw / 3.0

def collect_rotations(fchk_text, prog_options):
    mw = np.sum(_array_from_fchk(fchk_text, array_name="Real atomic weights"))
    rots = []
    if prog_options is not None:
        wls = prog_options.get('omega')
        if wls:
            # strip off the unit
            if isinstance(wls[-1], str):
                wls = wls[:-1]
            # better than doing one additional conversion just use the values in the fchk file
            au_freqs = _array_from_fchk(fchk_text, array_name="Frequencies for FD properties")
            # sort wls in acend. order, then swap so they are in aced energy order (same as freqs from fchk)
            wls = list(sorted(wls))
            wls.reverse()
            if (au_freqs is not None) and len(au_freqs) == len(wls):
                all_rot_tensors = _array_from_fchk(fchk_text, array_name="FD Optical Rotation Tensor").reshape(-1, 9)
                for w_au, w_nm, rot_tensor in zip(au_freqs, wls, all_rot_tensors):
                    val = compute_rotation(w_au, rot_tensor.reshape(3,3), mw)
                    rots.append({'value': val, 'wavelength': w_nm, 'gauge': 'GIAO'})
    return rots

def collect_hessian(fchk_text, natom):
    hessian_data = _array_from_fchk(fchk_text, array_name="Cartesian Force Constants")
    full_data = np.zeros((3*natom, 3*natom))
    test_data = np.zeros((3*natom, 3*natom))
    ut_index = np.triu_indices_from(full_data)
    lt_index = np.tril_indices_from(full_data)
    full_data[ut_index] = hessian_data
    full_data[lt_index] = hessian_data
    if np.allclose(test_data, full_data):
        return None
    else:
        return full_data

def collect_gradient(fchk_text):
    grad_data = _array_from_fchk(fchk_text, array_name="Cartesian Gradient")
    if np.allclose(grad_data, np.zeros_like(grad_data)):
        return None
    else:
        return grad_data

# execution steps

def exe_g09():
    sp.run('g09 input.com output.log', check=True, shell=True, env=os.environ)

def exe_fchk():
    sp.run(['formchk','-3','vices.chk','vices.fchk'], env=os.environ)

# input helpers
def _link_header(n):
    ret = []
    ret.append("--Link{}--".format(n))
    ret.append("%chk=vices")
    ret.append("%mem={}mb".format(cpu_info.memory()))
    ret.append("%nproc={}".format(cpu_info.ncore()))
    return ret

def _geometry_sec(input_json):
    mol_json = input_json.get('molecule')
    geom_array = mol_json.get('geometry')
    ret = []
    ret.append("{} {}".format(int(mol_json.get('charge')), mol_json.get('multiplicity')))
    for i_atom, symb in enumerate(mol_json.get('symbols')):
        ret.append("{} {:20.12f} {:20.12f} {:20.12f}".format(symb,
            geom_array[3*i_atom]   * physconst['bohr2angstroms'],
            geom_array[3*i_atom+1] * physconst['bohr2angstroms'],
            geom_array[3*i_atom+2] * physconst['bohr2angstroms']))
    ret.append('')
    return ret

def _title(driver, mol_name, modelchem_name):
    # blank line after title
    if mol_name is None:
        mol_name = "(no name)"
    if modelchem_name is None:
        modelchem_name = '(no name)'
    return ['{} mol {} mc {}'.format(driver,mol_name, modelchem_name), '']

def _rotation_wls(omegas):
    return ['{}nm'.format(x) for x in omegas]

def _end_input():
    return ["\n","\n", "\n"]

def _route(driver,input_json):
    route_line = ["#P"]
    if  driver == 'gradient':
        route_line.append("Force")
    elif driver == 'hessian':
        route_line.append("Freq")
    elif driver == 'rotation':
        route_line.append("")
    route_line.append("{}/{}".format(input_json['modelchem'].get('method'), input_json['modelchem'].get('basis')))
    route_line.append("scf(conver=11) int=ultrafine")
    if driver == 'rotation':
        route_line.append('polar=OptRot')
        route_line.append('cphf(RdFreq,conver=11)')

    # Blank line follows route section
    return [" ".join(route_line), '']


def write_input(driver, input_json):
    lines = []
    lines.extend(_link_header(0))
    lines.extend(_route(driver, input_json))
    lines.extend(_title(driver, input_json['molecule'].get('name'), input_json['modelchem'].get('name')))
    lines.extend(_geometry_sec(input_json))
    if driver == 'rotation':
        po = input_json.get('modelchem').get('program_options')
        if po:
            omegas = po.get('omega')
            if omegas:
                omegas = omegas[:-1]
                lines.extend(_rotation_wls(omegas))
    lines.extend(_end_input())
    infile_path = Path('input.com')
    infile_path.write_text("\n".join(lines))

# main driver
def run(input_json):
    calc_type = input_json['modelchem'].get('driver')
    output_json = {}
    output_json['raw_output'] = {}
    output_json['success'] = False
    output_json['output'] = {}

    try:
        # write the input file
        write_input(calc_type, input_json)
        # exe g09
        exe_g09()
        # exe fchk
        exe_fchk()

        # gather up stuff
        fchk_path = Path('vices.fchk')
        raw_out_text = Path('output.log').read_text()
        hess = collect_hessian(fchk_path.read_text(), len(input_json['molecule']['symbols']))
        grad = collect_gradient(fchk_path.read_text())
        rotations = collect_rotations(fchk_path.read_text(), input_json['modelchem'].get('program_options'))

        # store the raw output
        output_json['raw_output']['log'] = raw_out_text
        # the fchk file can be too large, mongo documents have a max size of 16mb so if the file is more than 12 we wont
        # include it, just tell where the path is
        if (fchk_path.stat().st_size // 1024**2) <= 12:
            output_json['raw_output']['fchk'] = fchk_path.read_text()
        else:
            output_json['raw_output']['fchk_path'] = str(fchk_path.resolve())
        if hess is not None:
            output_json['output']['hessian'] = util.pack_json_field(hess)
        if grad is not None:
            output_json['output']['gradient'] = util.pack_json_field(grad)
        if rotations:
            output_json['output']['rotations'] = rotations
        output_json['success'] = True
        return output_json
    except Exception as e:
        output_json['success'] = False
        output_json['raw_output']['error'] = "\n".join(traceback.format_exception(*sys.exc_info()))
        return output_json
