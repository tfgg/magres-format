#!python
import sys, os
import shutil
from magres.utils import find_all 

cwd = sys.argv[1]
target_dir = sys.argv[2]

file_paths = find_all(cwd, ".sh") + find_all(cwd, ".cell") + find_all(cwd, ".param") + find_all(cwd, ".DAT")

for path in file_paths:
  bits = path.split('/')

  for i in range(len(bits)):
    new_path = os.path.join(target_dir, *bits[:i])

    if not os.path.isdir(new_path):
      os.mkdir(new_path)

  new_path = os.path.join(target_dir, *bits)

  shutil.copy(path, new_path)

