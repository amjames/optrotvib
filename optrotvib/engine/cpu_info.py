import psutil
import cpuinfo
import socket
import os

_known_hn_sysn_coversions = {
        'nr': 'newriver',
        'br': 'blueridge',
        'dt': 'dragonstooth',
        'ca': 'cascades'
        }

def hostname():
    return socket.gethostname()

def cluster_name():
    cn = os.environ.get('SYSNAME')
    if cn:
        return cn
    cn = hostname()
    for pref, fulln in _known_hn_sysn_coversions.items():
        if cn.startswith(pref):
            return fulln
    return cn

def ncore():
    return psutil.cpu_count(logical=False)

def memory():
    "reasonable estimate of how much memory a job could use in MB"
    return int(psutil.virtual_memory().available * 0.9 / (1024*1024))

