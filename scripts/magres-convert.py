#!python
# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
import re
import math
import numpy
from magres.oldmagres import OldMagres

if __name__ == "__main__":
    magres_file = None
    castep_file = None

    if len(sys.argv) >= 2:
        try:
            magres_file = open(sys.argv[1]).read()
        except IOError:
            print("Could not load magres file '%s'" % sys.argv[1], file=sys.stderr)
            sys.exit(1)

    if len(sys.argv) >= 3:
        try:
            castep_file = open(sys.argv[2]).read()
        except IOError:
            print("Could not load castep output file '%s'" % sys.argv[2], file=sys.stderr)
            sys.exit(1)

    if castep_file is None:
        print(
            "WARNING: You have not supplied a .castep file as the second argument, so I don't know the crystal "
            "lattice.",
            file=sys.stderr)

    if magres_file is not None:
        oldmagres_file = OldMagres(magres_file, castep_file)

        magres_file = oldmagres_file.as_new_format()
        print(magres_file)
