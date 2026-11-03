"""Exact endpoint and inversion algebra through weight six.

This is calculation code, not a formal proof. Endpoint MZVs are reduced by
shuffle/stuffle relations, duality, and the depth-one zeta values. No fitted
integral result or selected Laurent coefficient enters this algebra.
"""
from collections import Counter
from functools import lru_cache
from fractions import Fraction as Q
import sympy as sp

Z2, Z3, Z5 = sp.symbols('Z2 Z3 Z5')

@lru_cache(None)
def shuffles(a, b):
    if not a:
        return {b: 1}
    if not b:
        return {a: 1}
    out = Counter()
    for w, c in shuffles(a[1:], b).items():
        out[(a[0],) + w] += c
    for w, c in shuffles(a, b[1:]).items():
        out[(b[0],) + w] += c
    return dict(out)

def word(indices):
    return tuple(a for n in indices for a in (0,) * (n - 1) + (1,))

def indices(w):
    out, n = [], 1
    for a in w:
        if a == 0:
            n += 1
        else:
            out.append(n)
            n = 1
    assert n == 1 and out and out[0] >= 2
    return tuple(out)

@lru_cache(None)
def compositions(n):
    if n == 0:
        return ((),)
    return tuple((k,) + rest for k in range(1, n + 1)
                 for rest in compositions(n - k))

@lru_cache(None)
def stuffle(a, b):
    if not a:
        return {b: 1}
    if not b:
        return {a: 1}
    out = Counter()
    for w, c in stuffle(a[1:], b).items():
        out[(a[0],) + w] += c
    for w, c in stuffle(a, b[1:]).items():
        out[(b[0],) + w] += c
    for w, c in stuffle(a[1:], b[1:]).items():
        out[(a[0] + b[0],) + w] += c
    return dict(out)

@lru_cache(None)
def regularized_word(w):
    """G(w;1), with G(0;1)=Reg G(1;1)=0; return MZV symbols."""
    if not w:
        return {(): Q(1)}
    if all(a == w[0] for a in w):
        return {}
    if w[-1] == 0:
        terms = shuffles(w[:-1], (0,))
    elif w[0] == 1:
        terms = shuffles((1,), w[1:])
    else:
        return {indices(w): Q((-1) ** w.count(1))}
    diagonal = terms[w]
    out = Counter()
    for other, c in terms.items():
        if other != w:
            for ix, value in regularized_word(other).items():
                out[ix] -= Q(c, diagonal) * value
    return {ix: c for ix, c in out.items() if c}

@lru_cache(None)
def mzv_values(n):
    if n == 0:
        return {(): sp.Integer(1)}
    assert 2 <= n <= 6
    admissible = tuple(ix for ix in compositions(n) if ix[0] >= 2)
    symbols = sp.symbols('v:' + str(len(admissible)))
    lookup = dict(zip(admissible, symbols))
    def endpoint(w):
        return sum(c * lookup[ix] for ix, c in regularized_word(w).items())
    depth_one = {2: Z2, 3: Z3, 4: 2 * Z2**2 / 5,
                 5: Z5, 6: 8 * Z2**3 / 35}
    equations = [lookup[(n,)] - depth_one[n]]
    for ix in admissible:
        # Iterated integral duality after u -> 1-u.
        dual = tuple(1 - a for a in reversed(word(ix)))
        equations.append(lookup[ix] - lookup[indices(dual)])
    for na in range(1, n):
        nb = n - na
        aa = ((1,),) if na == 1 else tuple(ix for ix in compositions(na) if ix[0] >= 2)
        bb = ((1,),) if nb == 1 else tuple(ix for ix in compositions(nb) if ix[0] >= 2)
        for a in aa:
            for b in bb:
                if a == (1,) and b == (1,):
                    continue
                sh = sum(c * (-1)**(len(a) + len(b)) * endpoint(w)
                         for w, c in shuffles(word(a), word(b)).items())
                st = sum(c * lookup.get(ix, sp.Symbol('div_' + str(ix)))
                         for ix, c in stuffle(a, b).items())
                # The divergent first index has the same coefficient and cancels.
                if a == (1,) or b == (1,):
                    # Shuffle regularization already set zeta(1)=0. The
                    # corresponding divergent stuffle coordinate equals the
                    # same shuffle-divergent coordinate; compute by direct
                    # word regularization rather than discarding it.
                    st = 0
                    for ix, c in stuffle(a, b).items():
                        if ix[0] >= 2:
                            st += c * lookup[ix]
                        else:
                            st += c * (-1)**len(ix) * endpoint(word(ix))
                    equations.append(sh - st)
                else:
                    product = mzv_values(na)[a] * mzv_values(nb)[b]
                    equations.extend((sh - product, st - product))
    solution = sp.linsolve(equations, symbols)
    assert len(solution) == 1, (n, solution)
    values = tuple(solution)[0]
    assert not any(set(v.free_symbols) & set(symbols) for v in values), (n, values)
    return dict(zip(admissible, map(sp.expand, values)))

@lru_cache(None)
def reginf_word(w):
    """Constant at +infinity: z=u/(1-u), with the real positive germ."""
    terms = {(): Q(1)}
    for a in w:
        assert a in (0, -1)
        choices = ((0, 1), (1, -1)) if a == 0 else ((1, -1),)
        nxt = Counter()
        for prefix, c in terms.items():
            for letter, sign in choices:
                nxt[prefix + (letter,)] += c * sign
        terms = nxt
    if not w:
        return sp.Integer(1)
    mzvs = mzv_values(len(w)) if len(w) > 1 else {}
    return sp.expand(sum(c * v * (1 if not ix else mzvs[ix])
                         for transformed, c in terms.items()
                         for ix, v in regularized_word(transformed).items()))

if __name__ == '__main__':
    for n in range(2, 7):
        print(n, len(mzv_values(n)))
