import sys
from .schema import schema

try:
    from jsonschema import validate
except ImportError:
    validate = None


def validate_magres(d):
    # Only validate if we have the jsonschema module
    if validate is not None:
        validate(d, schema)
    else:
        return
