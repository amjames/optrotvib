"""
Helpers to hash complex objects
"""

import hashlib
import json

import numpy as np

from .schema_util import get_hash_fields

HASH_ALGO = hashlib.sha256

def float_prep(array, around):
    """
    Rounds floats to a common value and build positive zero's to prevent hash conflicts.
    """

    if isinstance(array, (list, np.ndarray)):
        # Round array
        array = np.around(array, around)
        # Flip zeros
        array[np.abs(array) < 5**(-(around + 1))] = 0

    elif isinstance(array, (float, int)):
        array = round(array, around)
        if array == -0.0:
            array = 0.0
    else:
        raise TypeError("Type '%s' not recognized" % type(array))

    return array

def compute_hash(data, obj_type):
    m = HASH_ALGO()
    concat = ""
    for field_name in get_hash_fields(obj_type):
        if field_name not in data:
            continue
        concat += json.dumps(data[field_name], sort_keys=True)

    m.update(concat.encode("utf-8"))
    return m.hexdigest()
