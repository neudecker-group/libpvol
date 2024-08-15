from . import _xhcfflib as bindings

def initialize_calculator(*args, **kwargs):
    return bindings.Calculator(*args, **kwargs)

