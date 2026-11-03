"""Closed arithmetic reader and ordered-domain AST adapter; no evaluation of source code."""
import ast
import copy
from .color import require

REPRESENTATION = 'ordered-parameters/1'

def q(n): return dict(kind='rational', numerator=str(n), denominator='1')
def v(n): return dict(kind='var', name=n)
def op(n, *args): return dict(kind='op', op=n, args=list(args))

def arithmetic(text, names):
    def read(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return q(node.value)
        if isinstance(node, ast.Name) and node.id in names:
            return v(node.id)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return op('neg', read(node.operand))
        if isinstance(node, ast.BinOp):
            operations = {ast.Add:'add', ast.Sub:'sub', ast.Mult:'mul', ast.Div:'div', ast.Pow:'rpow'}
            require(type(node.op) in operations, 'unsupported arithmetic operator')
            return op(operations[type(node.op)], read(node.left), read(node.right))
        raise ValueError('unsupported or unbound arithmetic expression')
    return read(ast.parse(text.replace('^', '**'), mode='eval').body)

def integral_ast(j, rename_z=False):
    require(j.get('representation') == REPRESENTATION, 'ordered representation required')
    names = j['variables']; order = j['integration_order_outer_to_inner']
    require(len(names) == len(set(names)) and set(names) == set(order) and len(order) == len(names),
            'invalid ordered variable inventory')
    require(set(j['bounds']) == set(names), 'ordered bounds inventory mismatch')
    require(j.get('measure_semantics') == 'iterated_real_interval_integral', 'unsupported measure semantics')
    require(j.get('delta_constraints') == {}, 'unsupported delta constraints')
    factors = j['propagator_polynomials']
    require(list(factors) == [f'f{i}' for i in range(1,len(factors)+1)], 'propagator inventory mismatch')
    expected = f'({j["numerator"]}) * '+' * '.join(f'({f})^(-1+eps)' for f in factors.values())
    require(j['expanded'] == expected, 'expanded integrand differs from declared factors')
    require(j['integrand'] == f'({j["numerator"]}) * '+' * '.join(f'{f}^(-1+eps)' for f in factors),
            'integrand differs from declared factors')
    scope = {'z1', 'eps'}; bounds = {}
    for name in order:
        require(name not in scope, 'ordered binder shadows parameter')
        low, high = j['bounds'][name]
        bounds[name] = (arithmetic(low, scope), arithmetic(high, scope))
        scope.add(name)
    term = arithmetic(j['expanded'], scope)
    for name in reversed(order):
        low, high = bounds[name]
        term = dict(kind='interval_integral', variable=dict(name=name,type='Real'),
                    lower=low,upper=high,body=term)
    def rename(t):
        if isinstance(t, list): return [rename(x) for x in t]
        if not isinstance(t, dict): return t
        if t.get('kind') == 'var' and t['name'] == 'z1': return v('z')
        return {k:rename(x) for k,x in t.items()}
    return rename(term) if rename_z else copy.deepcopy(term)
