


def _make_job_json(molecule, job_spec, name):
    job_json = {}
    job_json['name'] = name
    job_json['modelchem'] = job_spec
    job_json['molecule'] = molecule.to_json()
    return job_json


def generate_jobs_for_stencil(stencil, job_spec, name):
    """make a directory for a job set, and create job jsons for each mode"""

    job_set_dir = Path() / name
    if job_set_dir.exists():
        return False
    job_set_dir.mkdir(parents=True)
    jobs = []
    jobs.append(_make_job_json(stencil['eq_molecule'], job_spec, "eq"))
    for mode_nm, mode_mol in stencil['modes'].items():
        jobs.append(_make_job_json(mode_mol, job_spec, mode_nm))

    return job_set_dir, jobs
