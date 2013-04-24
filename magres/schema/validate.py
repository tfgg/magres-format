import sys
from schema import schema
from jsonschema import validate

def validate_magres(d):
  validate(d, schema)

