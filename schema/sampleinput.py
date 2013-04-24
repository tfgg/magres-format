import json

input = {
    'atoms': {
      'lattice': [
        [[3.0, 0.0, 0.0],
         [0.0, 3.0, 0.0],
         [0.0, 0.0, 3.0]]
      ],
      'atom': [
        {
          "species": "Si",
          "label": "Si",
          "index": 1,
          "position": [0.0, 0.0, 0.0]
        },
        {
          "species": "Al",
          "label": "Al",
          "index": 15,
          "position": [1.5, 1.5, 1.5]
        }
      ]
    },
    'magres': {
      'efg_local': [
        {
          "atom": {
            "label": "Si",
            "index": 1
          }, 
          "V": [[0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]]
        }
       ],
      'isc': [
        {
          "atom1": {
            "label": "Si",
            "index": 1
          }, 
          "atom2": {
            "label": "Al",
            "index": 15
          }, 
          "K": [[0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]]
        }
      ]
    }
  }

if __name__ == "__main__":
  print json.dumps(input)
