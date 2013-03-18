#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import re
import math
import numpy

from magres.oldmagres import OldMagres

if __name__ == "__main__":
  if len(sys.argv) == 2:
    oldmagres_file = OldMagres(open(sys.argv[1]).read())

  if len(sys.argv) > 2:
    oldmagres_file = OldMagres(open(sys.argv[1]).read(), open(sys.argv[2]).read())

  magres_file = oldmagres_file.as_new_format()

  print magres_file
