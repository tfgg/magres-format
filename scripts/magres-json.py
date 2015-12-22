#!python
from __future__ import print_function
import sys
import magres.format

if __name__ == "__main__":
    magres_file = magres.format.MagresFile(open(sys.argv[1]))

    print(magres_file.as_json())
