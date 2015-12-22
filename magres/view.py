from numpy import mean

def is_iter(xs):
  if type(xs) in [ListPropertyView, list]:
    return True
  else:
    return False

def flatten(xs):
  if all(map(is_iter,xs)):
    if len(xs) > 0:
      xs = list(map(flatten, xs))
      return sum(xs[1:], xs[0])
    else:
      return ListPropertyView([])
  else:
    return xs

class ListPropertyView(list):
  """
    Allows property accessors on lists of objects. E.g.

    x = [A(1), A(2), A(3)]
    x.a = [1,2,3]

    if A(a) is an object that stores a in self.a
  """

  def mean(self, *args, **kwargs):
    return mean([x for x in self], *args, **kwargs) 

  def __getattr__(self, prop):
    if any([hasattr(x, prop) for x in self]):
      return ListPropertyView(flatten([getattr(x, prop) for x in self if hasattr(x, prop)]))
    else:
      #raise AttributeError("Property '{}' not present".format(prop))
      return ListPropertyView([])

  def __repr__(self):
      return "ListPropertyView([{}])".format(", ".join(repr(x) for x in self))

  #def _repr_html_(self):
  #  return html_repr.list_view(self)

  def __call__(self, *args, **kwargs):
    return flatten(ListPropertyView([x(*args, **kwargs) for x in self]))

  def __add__(self, b):
    return ListPropertyView(super(ListPropertyView, self).__add__(b))


