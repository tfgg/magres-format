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
