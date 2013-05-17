class lazyproperty(object):
  '''
    Meant to be used for lazy evaluation of an object attribute.
    Property should represent non-mutable data, as it replaces itself.
  '''

  def __init__(self, fget=None, fset=None, fdel=None, doc=None):
    self.fget = fget
    self.func_name = fget.__name__

    self.__doc__ = doc

    if doc is None and fget is not None:
      self.__doc__ = fget.__doc__

    self.__name__ = fget.__name__

  def __get__(self,obj,cls):
    # When the property on the owner instance is first accessed we 
    # calculate its value and self-modify to insert the value on 
    # the owner instance
    
    if obj is None:
      return self

    return self.fget(obj)

