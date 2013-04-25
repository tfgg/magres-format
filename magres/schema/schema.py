import json
import sys
from jsonschema import validate

schema = {
            "title": "magres file",
            "description": "JSON serialisation of a ab-intio magres file",
            "type": "object",
            "properties": {
              "calculation": {
                "type": "object",
                "patternProperties": {
                  ".*": {
                    "type": "array",
                    "items": {
                      "type": "array",
                      "items": {
                        "type": "string"
                      }
                    }
                  }
                }
              },
              "atoms": {
                "type": "object",
                "properties": {
                  "lattice": {
                    "type": "array",
                    "items": {"$ref": "#/definitions/tensor33"},
                    "maxItems": 1,
                    "minItems": 0
                  },
                  "atom": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "species": {
                          "type": "string"
                        },
                        "label": {
                          "type": "string"
                        },
                        "index": {
                          "type": "integer"
                        },
                        "position": {
                          "type": "array",
                          "items": {"type": "number"},
                          "maxItems": 3,
                          "minItems": 3
                        }
                      }
                    }
                  }
                }
              },
              "magres": {
                "type": "object",
                "patternProperties": {
                  "(efg|efg_local|efg_ions|efg_nonlocal)": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "atom": {
                          "$ref": "#/definitions/atomindex"
                        },
                        "V": {
                          "$ref": "#/definitions/tensor33"
                        }
                      }
                    }
                  },
                  "(isc|isc_spin|isc_fc|isc_orbital_p|isc_orbital_d)": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "atom1": {
                          "$ref": "#/definitions/atomindex"
                        },
                        "atom2": {
                          "$ref": "#/definitions/atomindex"
                        },
                        "K": {
                          "$ref": "#/definitions/tensor33"
                        }
                      }
                    }
                  }
                }
              }
            },
            "required": ["magres"],
            "definitions": {
              "atomindex": {
                "type": "object",
                "properties": {
                  "label": {
                    "type": "string",
                  },
                  "index": {
                    "type": "integer",
                    "minimum": 0
                  }
                }
              },
              "tensor33": {
                "type": "array",
                "items": {"type": "array", "items": {"type": "number"}, "maxItems": 3, "minItems": 3},
                "maxItems": 3,
                "minItems": 3
              }
            }
         }

#input = json.load(open(sys.argv[1]))

#validate(input, schema)

if __name__ == "__main__":
  print json.dumps(schema)

