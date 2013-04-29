import os

def insideout():
  """
    Count up in positive numbers and down in negative numbers
  """

  yield 0

  i = 1

  while True:
    yield i
    yield -i
    i += 1

def find_all_magres(dir):
  calcs = []
  for f in os.listdir(dir):
    path = os.path.join(dir, f)
    if ".magres" in f:
      calcs.append(path)
    elif os.path.isdir(path):
      calcs += find_all_magres(path)
  
  return calcs
