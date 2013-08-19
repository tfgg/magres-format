import sys, os
from magres.utils import find_all

files = find_all('orig', '.magres') + find_all('low', '.magres')

for file in files:
  lines = []

  new_path = os.path.join("processed", file)
  
  new_file = open(new_path, "w+")

  for line in open(file):
    cols = line.split()

    if len(cols) < 1:
      continue

    if cols[0] == "atom":
      cols[2] = cols[1]
      print >>new_file, " ".join(cols)
    else:
      print >>new_file, line.strip()
