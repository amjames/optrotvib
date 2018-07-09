import os
import subprocess as sp
import psutil
import time
from lxml import etree
from pathlib import Path

class temporary_move(object):
    def __init__(self):
        self.return_dir = None

    def __enter__(self, t_dir):
        self.return_dir = Path.cwd()
        os.chdir(str(t_dir))

    def __exit__(self, *exc):
        os.chdir(str(self.return_dir))


def qsub(args):
    full_args = "-W group_list={} {}".format(os.environ.get("SYSNAME"), args)
    proc = sp.run("qsub {} submit.pbs".format(full_args), shell=True, check=True, universal_newlines=True, stdout=sp.PIPE)
    return proc.stdout.split('.')[0]

def qstat(job_id):
    proc = sp.run("qstat -x {}".format(job_id), shell=True, stdout=sp.PIPE, universal_newlines=True, check=False)
    # if the job is not found then we will assume it was completed and now is missing
    if proc.returncode != 0:
        return 'C'
    state = etree.XML(proc.stdout.strip()).find("Job").find('job_state').text
    logger.info("QSTAT: job {} STATE = '{}'".format(job_id, state))
    return state

def write_submit_script(env, job_path):
    lines = env + ['optrotvib.engine']
    Path(job_path / 'submit.pbs').write_text("\n".join(lines))

def submit(compute_env, qsub_args, compute_dir, job):
    job_path = compute_dir / job['name']
    if not job_path.exists():
        job_path.mkdir(parents=True)

    write_submit_script(env, job_path)
    Path(job_path / 'input.json').write_text(json.dumps(job))
    with temporary_move(job_path):
        jid = qsub(qsub_args)
        Path('job.id').write_text(jid)
    return jid

def check_by_dir(job_dir):
    jid = Path(job_dir / 'job.id').read_text()
    return qstat(jid)

def check_by_id(jid):
    return qstat(jid)
