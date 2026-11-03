"""Exact fixed-basis color algebra; depends only on Python and SymPy."""
import ast
import sympy as sp
A,F=sp.symbols('CA CF')
REL=A*A-1-2*A*F
class ColorError(ValueError):pass

def require(condition,message):
    if not condition:raise ColorError(message)

def q(n):return dict(kind='rational',numerator=str(n),denominator='1')
def v(name):return dict(kind='var',name=name)
def op(name,*args):return dict(kind='op',op=name,args=list(args))

def parse(text):
    require(isinstance(text,str) and len(text)<=2048,'color expression length')
    try:tree=ast.parse(text.replace('^','**'),mode='eval')
    except SyntaxError as exc:raise ColorError('invalid color arithmetic') from exc
    require(sum(1 for _ in ast.walk(tree))<=256,'color expression budget')
    def degree(node):
        if isinstance(node,ast.Name):return 1
        if isinstance(node,ast.Constant):return 0
        if isinstance(node,ast.UnaryOp):return degree(node.operand)
        require(isinstance(node,ast.BinOp),'only arithmetic expressions are supported')
        if isinstance(node.op,ast.Pow):
            try:exponent=ast.literal_eval(node.right)
            except (ValueError,TypeError):raise ColorError('literal integer exponent required')
            require(type(exponent) is int and abs(exponent)<=8,'exponent budget')
            return degree(node.left)*abs(exponent)
        left,right=degree(node.left),degree(node.right)
        return max(left,right) if isinstance(node.op,(ast.Add,ast.Sub)) else left+right
    require(degree(tree.body)<=16,'color polynomial degree budget')
    def visit(node):
        if isinstance(node,ast.Name):
            require(node.id in ['CA','CF'],'only CA and CF are permitted')
            return (A if node.id=='CA' else F),v(node.id)
        if isinstance(node,ast.Constant):
            require(type(node.value) is int and abs(node.value)<10**30,'exact integer literals required')
            return sp.Integer(node.value),q(node.value)
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,(ast.USub,ast.UAdd)):
            x,t=visit(node.operand)
            return (-x,op('neg',t)) if isinstance(node.op,ast.USub) else (x,t)
        require(isinstance(node,ast.BinOp),'only arithmetic expressions are supported')
        x,tx=visit(node.left);y,ty=visit(node.right)
        if isinstance(node.op,ast.Add):return x+y,op('add',tx,ty)
        if isinstance(node.op,ast.Sub):return x-y,op('sub',tx,ty)
        if isinstance(node.op,ast.Mult):return x*y,op('mul',tx,ty)
        if isinstance(node.op,ast.Div):
            den=sp.cancel(y)
            coeff,power=den.as_coeff_exponent(A)
            require(coeff.is_Rational and coeff!=0 and power.is_Integer and power>=0,'denominators must be nonzero rational multiples of CA powers')
            return x/y,op('div',tx,ty)
        require(isinstance(node.op,ast.Pow) and y.is_Integer and -8<=y<=8,'only integer color powers -8..8')
        if y<0:
            require(x==A,'negative powers only for CA')
        t=q(1)
        for _ in range(abs(int(y))):t=op('mul',t,tx)
        return x**y,(op('div',q(1),t) if y<0 else t)
    return visit(tree.body)


def project(n,text):
    require(type(n) is int and 1<=n<=8,'loop order must be 1..8')
    expr,term=parse(text)
    basis=[A**(n-k)*F**k for k in range(1,n+1)]
    unknown=sp.symbols('a0:'+str(n))
    reduced=sp.cancel((expr-sum(c*b for c,b in zip(unknown,basis))).subs(F,(A*A-1)/(2*A)))
    polynomial=sp.Poly(sp.fraction(reduced)[0],A)
    matrix,rhs=sp.linear_eq_to_matrix(polynomial.all_coeffs(),unknown)
    require(matrix.rank()==n,'basis coordinates are not uniquely determined')
    solutions=sp.linsolve((matrix,rhs),unknown)
    require(solutions is not sp.EmptySet and len(solutions)==1,'expression is outside the fixed-loop basis')
    coordinates=list(next(iter(solutions)))
    require(all(c.is_Rational for c in coordinates),'color coordinates must be rational constants')
    reconstruction=sp.expand(sum(c*b for c,b in zip(coordinates,basis)))
    numerator,denominator=sp.fraction(sp.cancel(expr-reconstruction))
    quotients,remainder=sp.reduced(numerator,[REL],A,F)
    require(remainder==0,'color reconstruction has a nonzero residual')
    return dict(expression=term,coordinates=coordinates,reconstruction=reconstruction,
                quotient=quotients[0] if quotients else sp.Integer(0),denominator=denominator,basis=basis)
