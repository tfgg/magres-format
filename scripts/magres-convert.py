#!python
# -*- coding: utf-8 -*-
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
      print >>sys.stderr, "Could not load magres file '%s'" % sys.argv[1]
      sys.exit(1)

  if len(sys.argv) >= 3:
    try:
      castep_file = open(sys.argv[2]).read()
    except IOError:
      print >>sys.stderr, "Could not load castep output file '%s'" % sys.argv[2]
      sys.exit(1)

  if magres_file is not None:
    oldmagres_file = OldMagres(magres_file, castep_file)

    magres_file = oldmagres_file.as_new_format()
    print magres_file
