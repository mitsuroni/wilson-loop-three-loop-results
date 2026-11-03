"""Single-use prepared C to Wilson I normalization and scalar Laurent output."""
from fractions import Fraction
import sympy as sp
from .algebra import rational, from_sympy, gamma_series


def prepared_normalization(data):
    if data['kind']!='prepared': raise ValueError('explicit I assembly is not yet supported')
    orbit=data['orbit']; n=data['loop_order']
    color=orbit['color_preparation']; euler=orbit['euler_preparation']
    a1=Fraction(color['a1']); sym=orbit['c4_symmetry']['orbit_multiplicity']
    if type(sym) is not int or sym<=0 or color['coordinates'][0]!=str(a1):raise ValueError('invalid color/symmetry coordinates')
    factor='Exp[-EulerGamma*eps]' if n==1 else f'Exp[-{n}*EulerGamma*eps]'
    expected=(f"({color['projected_color']}) * ({factor}) * ((4*Pi*alphaS)^{n})"
              f" * (Gamma[1-eps]^{n}/(4*Pi^(2-eps))^{n}) * ((-2*v14+I0)^({n}*eps)) * ((-2)^(-{n}))")
    if orbit['C']['combined']!=expected or euler['gamma_power']!=n or euler['factor']!=factor:
        raise ValueError('prepared C does not match the single-Euler assembly recipe')
    if orbit['C']['components'].get('euler_compensation')!=factor:raise ValueError('Euler metadata mismatch')
    return {'format':'wilson-2ne-assembly/0.1','n':n,'kappa':str(Fraction(sym)*a1/(2*(-2)**n)),
            'euler':'consumed from prepared C exactly once',
            'formula':'kappa*G(eps)^n*(J(z)+z^(n*eps)*J(1/z))/(1+z^(n*eps))'}


def physical_expression(j,normalization,eps,z):
    n=normalization['n']; k=sp.Rational(normalization['kappa'])
    g=sp.exp(-sp.EulerGamma*eps)*sp.gamma(1-eps)
    return k*g**n*(j+z**(n*eps)*j.subs(z,1/z))/(1+z**(n*eps))


def scalar_hlog(expr,z):
    """Encode polynomials in log(z) as actual zero-letter GPL words."""
    ell=sp.Symbol('_logz'); polynomial=sp.expand(expr).subs(sp.log(z),ell)
    try: poly=sp.Poly(polynomial,ell)
    except sp.PolynomialError as exc: raise ValueError('scalar expression outside log-polynomial rule') from exc
    terms=[]
    for (degree,),coefficient in poly.terms():
        if coefficient==0:continue
        coefficient=coefficient.xreplace({z:sp.Symbol('z')})
        terms.append({'coefficient':from_sympy(sp.expand(coefficient*sp.factorial(degree))),
                      'gpl':{'word':[rational(0)]*degree,'endpoint':{'kind':'var','name':'z'},
                             'basepoint':rational(0),'branch':'positive-real','endpoint_kind':'ordinary'},'periods':[]})
    return {'format':'wilson-hlog-sum/0.1','kernel':'dt/(t-a)','parameters':['z'],'terms':terms}


def expand_scalar(j,normalization,quantity,through,eps,z,loop_order):
    expression=j if quantity=='J' else physical_expression(j,normalization,eps,z)
    # Each adjacent propagator has at most two endpoint poles in the supported
    # elementary chain rule. Its concrete Gamma expression supplies the series.
    lower=-2*loop_order
    if through<lower:return {},lower
    jet=gamma_series(expression,eps,through-lower,output_through=through)
    coefficients={}
    for power in range(lower,through+1):
        coefficients[str(power)]=scalar_hlog(sp.expand(jet).coeff(eps,power),z)
    if any(t.has(eps) for t in (sp.expand(jet).coeff(eps,k) for k in range(lower,through+1))):
        raise ValueError('unexpanded regulator in a purported coefficient')
    return coefficients,lower
