import json
import sys
from jsonschema import validate

schema = {
            "title": "magres file",
            "description": "JSON serialisation of a ab-initio magres file",
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
                  "units": {
                    "$ref": "#/definitions/units"
                  },
                  "symmetry": {
                    "$ref": "#/definitions/symmetry"
                  },
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
                "additionalProperties": {
                  "type": "array",
                  "items": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  }
                },
                "patternProperties": {
                  "units": {
                    "$ref": "#/definitions/units"
                  },
                  "(efg|efg_.*)": {
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
                  "(isc|isc_.*)": {
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
                  },
                  "(ms|ms_.*)": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "atom": {
                          "$ref": "#/definitions/atomindex"
                        },
                        "sigma": {
                          "$ref": "#/definitions/tensor33"
                        }
                      }
                    }
                  }
                }
              }
            },
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
              },
              "units": {
                "type": "array",
                "items": {
                  "type": "array",
                  "items": {"type": "string"},
                  "maxItems": 2,
                  "minItems": 2
                }
              },
              "symmetry": {
                "type": "array",
                "items": {
                  "type": "string",
                }
              }
            }
         }

#input = json.load(open(sys.argv[1]))

#validate(input, schema)

if __name__ == "__main__":
  print json.dumps(schema)

