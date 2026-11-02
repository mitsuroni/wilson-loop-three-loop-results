"""Independent computational checks; none of these are Lean receipts."""
import functools
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.stats import qmc
from empty_math import eps, z, ell, SUBS, to_sympy, word_integral


def decimal(x):
    return mp.nstr(x, 50)


@functools.lru_cache(maxsize=None)
def check_word(word):
    checks = []
    with mp.workdps(70):
        for e in (mp.mpf(1)/2, mp.mpf(2)/3, mp.mpf(3)/4):
            def f(c,x):
                return x**(e-1) if c == "A" else (1-x)**(e-1)
            def left(c,x):
                return x**e/e if c == "A" else -mp.expm1(e*mp.log1p(-x))/e
            def right(c,x):
                return -mp.expm1(e*mp.log(x))/e if c == "A" else (1-x)**e/e
            def integrand(x):
                if x == 0 or x == 1: return mp.mpf(0)
                if len(word) == 1: return f(word[0],x)
                if len(word) == 2: return left(word[0],x)*f(word[1],x)
                assert len(word) == 3
                return left(word[0],x)*f(word[1],x)*right(word[2],x)
            observed = mp.quad(integrand,[0,mp.mpf("0.5"),1])
            exact = sp.lambdify(eps,word_integral(word).subs(SUBS),"mpmath")(e)
            error = abs(observed-exact)
            relative = error/max(abs(exact),mp.mpf(1))
            assert relative < mp.mpf("1e-25"), (word,e,error)
            checks.append({"eps":decimal(e),"quadrature":decimal(observed),
                           "closed_form":decimal(exact),"absolute_error":decimal(error),
                           "scaled_error":decimal(relative),"working_dps":70})
    return {"word":word,"method":"direct one-dimensional integral of the ordered word after elementary endpoint integration",
            "checks":checks,"passed":True}


def original_polynomial_checks(row, a):
    j = row["orbit"]["J"]
    original = to_sympy(j["expanded"],j["variables"])
    checks=[]
    for value in (1,2,3):
        expr = sp.expand(original.subs({eps:value,z:sp.Rational(2,3)}))
        steps=[]
        for name in reversed(j["integration_order_outer_to_inner"]):
            v = sp.Symbol(name,positive=True)
            low,high = [to_sympy(s,j["variables"]) for s in j["bounds"][name]]
            poly=sp.Poly(expr,v)
            expr=sp.expand(sum(c*(high**(power[0]+1)-low**(power[0]+1))/(power[0]+1)
                               for power,c in poly.terms()))
            steps.append({"integrated_variable":name,"remaining_expression":str(expr)})
        candidate = a["J_exact"].subs({eps:value,z:sp.Rational(2,3)})
        residual=sp.simplify(expr-candidate)
        assert residual == 0, (row["orbit_id"],value,residual)
        checks.append({"eps":value,"z":"2/3","original_iterated_value":str(expr),
                       "closed_value":str(candidate),"residual":str(residual),"steps":steps})
    return {"method":"exact polynomial integration of prepared expanded J in its original outer-to-inner order (integrated inside out)",
            "scope":"J only; epsilon=1,2,3 are not evaluations of the Gamma-prefactored I",
            "checks":checks,"passed":True}


def raw_domain_qmc(row,a):
    j = row["orbit"]["J"]
    names = j["integration_order_outer_to_inner"]
    checks=[]
    power=16
    repetitions=8
    for point,e in enumerate((0.5,0.7)):
        estimates=[]
        for rep in range(repetitions):
            seed=4100+point*100+rep
            u=qmc.Sobol(d=len(names),scramble=True,seed=seed).random_base2(power)
            samples={}
            for i,v in enumerate(names):
                x=u[:,i]**(1/e)
                samples[v]=x if a["letters"][v] == "A" else 1-x
            inside=np.ones(len(u),dtype=bool)
            for v,(low,high) in j["bounds"].items():
                assert low=="0"
                if high!="1": inside &= samples[v] <= samples[high]
            estimates.append(float(np.mean(inside)))
        mean=float(np.mean(estimates))
        se=float(np.std(estimates,ddof=1)/np.sqrt(repetitions))
        with mp.workdps(40):
            exact=float(sp.lambdify(eps,(eps**(2*row["loop_order"])*a["K"]).subs(SUBS),"mpmath")(mp.mpf(str(e))))
        delta=abs(mean-exact)
        tolerance=max(8*se,8/(2**power))
        assert delta <= tolerance,(row["orbit_id"],e,delta,tolerance)
        checks.append({"eps":e,"seeds":[4100+point*100+r for r in range(repetitions)],
                       "samples_per_scramble":2**power,"scrambles":repetitions,
                       "probability_estimates":estimates,"mean":mean,"standard_error_estimate":se,
                       "exact_probability":exact,"absolute_difference":delta,"acceptance_tolerance":tolerance})
    return {"method":"randomized Sobol expectation of the ORIGINAL bounds indicator under independent normalized A/B endpoint densities; J=z^(r eps)*eps^(-2n)*probability",
            "scope":"statistical consistency, not a certified numerical enclosure; acceptance max(8 SE,8/N)",
            "checks":checks,"passed":True}


def cauchy_coefficients(row,a,coefficients):
    records=[]
    n=row["loop_order"]
    with mp.workdps(70):
        fun=sp.lambdify((eps,ell),a["I_compact"].subs(SUBS),"mpmath")
        radius=mp.mpf("0.025")
        size=64
        roots=[mp.exp(2j*mp.pi*j/size) for j in range(size)]
        for zz in (mp.mpf("0.5"),mp.mpf(1),mp.mpf(2)):
            logz=mp.log(zz)
            values=[(radius*w)**(2*n)*fun(radius*w,logz) for w in roots]
            comparisons=[]
            for exponent in range(-2*n,1):
                k=exponent+2*n
                observed=sum(values[j]/roots[j]**k for j in range(size))/size/radius**k
                expected=sp.lambdify(ell,coefficients[exponent],"mpmath")(logz)
                error=abs(observed-expected)
                assert error/max(1,abs(expected)) < mp.mpf("1e-35"),(row["orbit_id"],exponent,error)
                comparisons.append({"power":exponent,"extracted_real":decimal(mp.re(observed)),
                                    "imaginary_residual":decimal(mp.im(observed)),
                                    "expected":decimal(expected),"absolute_error":decimal(error)})
            records.append({"z":decimal(zz),"coefficients":comparisons})
    return {"method":"Cauchy/Fourier coefficient extraction from the exact Gamma expression after multiplying by eps^(2n); independent of symbolic series recurrence",
            "working_dps":70,"circle_radius":"0.025","samples":64,"checks":records,"passed":True}


def remainder_samples(row,a,coefficients):
    records=[]
    with mp.workdps(80):
        fun=sp.lambdify((eps,ell),a["I_compact"].subs(SUBS),"mpmath")
        partial=sp.lambdify((eps,ell),sum(v*eps**k for k,v in coefficients.items() if k<=0),"mpmath")
        c1=sp.lambdify(ell,coefficients[1],"mpmath")
        for zz in (mp.mpf("0.5"),mp.mpf(1),mp.mpf(2)):
            logz=mp.log(zz)
            expected=c1(logz)
            previous=None
            points=[]
            for denominator in (1000,2000,4000):
                e=mp.mpf(1)/denominator
                observed=(fun(e,logz)-partial(e,logz))/e
                error=abs(observed-expected)
                if previous is not None:
                    assert error <= previous*mp.mpf("0.75")+mp.mpf("1e-40")
                previous=error
                points.append({"eps":str(e),"remainder_over_eps":decimal(observed),
                               "next_coefficient":decimal(expected),"difference":decimal(error)})
            records.append({"z":decimal(zz),"points":points})
    return {"method":"finite part remainder divided by eps converges to symbolic next coefficient; O(eps) itself follows from analytic expression",
            "note":"Cauchy check above covers through power 0; next coefficient here is symbolic, with numerical convergence only.",
            "working_dps":80,"checks":records,"passed":True}
