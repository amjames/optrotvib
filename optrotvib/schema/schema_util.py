"""
Assists in grabbing the requisite schema
"""

import copy
import glob
import os
import jsonschema
import json

from pathlib import Path
__all__ = ["get_schema", "validate", "get_hash_fields", "get_valid_fields", "get_indices"]

_schemas = {}
schemas_root = Path(__file__).parent
for schema_p in list(schemas_root.glob("*.schema.json")):
    schema_n = schema_p.name.split('.')[0]
    _schemas[schema_n] = json.loads(schema_p.read_text())

_definitions = json.loads(Path(schemas_root / 'definitions.json').read_text())
for sn in _schemas.keys():
    _schemas[sn]["definitions"] = copy.deepcopy(_definitions)


def get_valid_fields(name):
    if name not in _schemas.keys():
        raise KeyError("Schema name %s not found." % name)
    return list(_schemas[name]['properties'].keys())

def get_index(name):
    if name not in _schemas.keys():
        raise KeyError("Schema name %s not found." % name)
    return copy.deepcopy(_schemas[name].get("index", []))

def get_hash_fields(name):
    if name not in _schemas.keys():
        raise KeyError("Schema name %s not found." % name)
    return copy.deepcopy(_schemas[name]["hash_fields"])

def get_schema(name):
    if name not in _schemas.keys():
        raise KeyError("Schema name %s not found." % name)
    return copy.deepcopy(_schemas[name])

def validate(data, schema_name, return_errors=False):
    schema = get_schema(schema_name)
    error_gen = jsonschema.Draft4Validator(schema).iter_errors(data)
    errors = [x for x in error_gen]
    if len(errors):
        if return_errors:
            return errors
        else:
            error_msg = "Error validating schema '%s'!\n" % schema_name
            error_msg += "Data: \n" + json.dumps(data, indent=2)
            error_msg += "\n\nJSON Schema errors as follow:\n"
            error_msg += "\r".join(x.message for x in errors)
            error_msg += "\n"

            raise ValueError(error_msg)
    else:
        return True
