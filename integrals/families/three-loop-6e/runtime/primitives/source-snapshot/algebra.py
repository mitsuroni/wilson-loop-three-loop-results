"""Controlled conversion between the existing scalar AST and exact SymPy objects."""
from fractions import Fraction
import sympy as sp


def rational(value):
    value = Fraction(value)
    return {'kind': 'rational', 'numerator': str(value.numerator), 'denominator': str(value.denominator)}


def ast_to_sympy(node, names=None):
    if not isinstance(node, dict):
        raise ValueError('scalar AST required')
    kind = node.get('kind')
    if kind == 'rational' and set(node) == {'kind', 'numerator', 'denominator'}:
        n, d = int(node['numerator']), int(node['denominator'])
        if d <= 0 or str(n) != node['numerator'] or str(d) != node['denominator'] or Fraction(n, d).denominator != d:
            raise ValueError('noncanonical rational')
        return sp.Rational(n, d)
    if kind == 'var' and set(node) == {'kind', 'name'}:
        name = node['name']
        if not isinstance(name, str) or not name.isidentifier() or (names is not None and name not in names):
            raise ValueError('unbound scalar variable')
        return names[name] if isinstance(names, dict) else sp.Symbol(name)
    if kind != 'op' or set(node) != {'kind', 'op', 'args'} or not isinstance(node['args'], list):
        raise ValueError('unsupported scalar AST')
    functions = {'add': (2, lambda a,b:a+b), 'sub': (2, lambda a,b:a-b),
                 'mul': (2, lambda a,b:a*b), 'div': (2, lambda a,b:a/b),
                 'rpow': (2, lambda a,b:a**b), 'neg': (1, lambda a:-a),
                 'exp': (1,sp.exp), 'log': (1,sp.log), 'gamma': (1,sp.gamma),
                 'zeta': (1,sp.zeta), 'pi': (0,lambda:sp.pi),
                 'euler_gamma': (0,lambda:sp.EulerGamma)}
    if node['op'] not in functions:
        raise ValueError('unsupported scalar operation')
    arity, function = functions[node['op']]
    if len(node['args']) != arity:
        raise ValueError('scalar arity mismatch')
    return function(*(ast_to_sympy(a, names) for a in node['args']))


def from_sympy(expr):
    expr = sp.sympify(expr) if isinstance(expr, (int,Fraction)) else expr
    if expr.is_Rational:
        return rational(Fraction(int(expr.p), int(expr.q)))
    if expr == sp.pi:
        return {'kind':'op','op':'pi','args':[]}
    if expr == sp.EulerGamma:
        return {'kind':'op','op':'euler_gamma','args':[]}
    if expr.is_Symbol:
        return {'kind':'var','name':str(expr)}
    operations = {sp.Add:'add',sp.Mul:'mul',sp.Pow:'rpow',sp.exp:'exp',sp.log:'log',sp.gamma:'gamma',sp.zeta:'zeta'}
    if expr.func not in operations:
        raise ValueError('unsupported symbolic output: '+str(expr.func))
    op = operations[expr.func]
    args = [from_sympy(a) for a in expr.args]
    if op in ('add','mul'):
        result = args[0]
        for a in args[1:]: result = {'kind':'op','op':op,'args':[result,a]}
        return result
    return {'kind':'op','op':op,'args':args}


def freeze_math(value):
    if isinstance(value,sp.Basic): return from_sympy(value)
    if isinstance(value,Fraction): return rational(value)
    if isinstance(value,dict): return {str(k):freeze_math(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [freeze_math(v) for v in value]
    if isinstance(value,(str,int,float,bool)) or value is None:return value
    raise ValueError('unsupported calculation metadata type: '+type(value).__name__)


def gamma_series(expr, eps, through, *, output_through=None):
    """Exact finite Laurent expansion; pole budgets are propagated internally.

    ``output_through`` preserves the explicit output/working-order separation
    of existing callers. No general CAS series search is used here.
    """
    from .jets import finite_laurent
    target = through if output_through is None else output_through
    values = finite_laurent(expr, eps, target)
    return sp.Add(*(coefficient*eps**power for power,coefficient in values.items())).expand()
