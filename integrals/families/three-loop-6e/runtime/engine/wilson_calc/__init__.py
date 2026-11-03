"""Deterministic ordered Wilson-integral calculations, with explicit coverage."""
__version__='0.2.0'


def solve(request,*,store=None):
    from .solver import solve as _solve
    return _solve(request,store=store)


__all__=['solve','__version__']
