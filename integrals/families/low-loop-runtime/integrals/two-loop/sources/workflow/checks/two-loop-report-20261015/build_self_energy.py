"""Build the two-loop pre-read's self-energy coefficients through eps**2.

The current r002 full-epsilon results and established M-moment source are input.
Stored finite-order coefficients are comparison data, never construction data.
This is algebraic expansion/report generation, not a new integral calculation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
TEX = ROOT / "benchmarks/renormalization/two-loop/report/self-energy-results.tex"
e, z, ell = s.symbols("eps z ell")
M01, M11, M21, M03, M31, M13 = s.symbols("M01 M11 M21 M03 M31 M13")
LOCALS = {str(x): x for x in (e, z, ell, M01, M11, M21, M03, M31, M13)}
bindings: dict[str, str] = {}
checks: list[dict] = []


def bind(rel):
    path = ROOT / rel
    bindings[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read(rel):
    return json.loads(bind(rel).read_text(encoding="utf-8-sig"))


def check(name, lhs, rhs=0, kind="exact"):
    residual = s.factor(s.expand(lhs - rhs))
    good = residual == 0
    checks.append(dict(name=name, kind=kind, passed=good, residual=str(residual)))
    assert good, (name, residual)


def truth(name, good, kind="hash"):
    checks.append(dict(name=name, kind=kind, passed=bool(good)))
    assert good, name


def exp_series(logpoly, order):
    a = s.Poly(logpoly, e)
    b = [s.S.One]
    for n in range(1, order + 1):
        b.append(s.expand(sum(k * a.nth(k) * b[n-k] for k in range(1, n+1)) / n))
    return sum(b[n] * e**n for n in range(order + 1))


def canonical(x):
    return s.expand(x.subs(M03, -M01**2).subs(ell**2, -2*M01-s.pi**2))


def tex_expr(x):
    names = {ell: r"\ell", M01: r"M_{0,1}", M11: r"M_{1,1}",
             M21: r"M_{2,1}", M03: r"M_{0,3}", M31: r"M_{3,1}", M13: r"M_{1,3}"}
    out = s.latex(s.expand(x), symbol_names=names)
    for n in (3, 5):
        out = out.replace(r"\zeta\left("+str(n)+r"\right)", r"\zeta_{"+str(n)+"}")
    return out


def equation_lines(lhs, x, count=4):
    terms = s.expand(x).as_ordered_terms()
    parts = []
    for i in range(0, len(terms), count):
        text = tex_expr(sum(terms[i:i+count]))
        if i and not text.startswith("-"):
            text = "+ " + text
        parts.append((lhs + " &=" if i == 0 else " &") + text)
    return r"\\*[2pt]".join(parts)


def coeff_block(symbol, coefficients):
    lines = [equation_lines(symbol+"_{"+str(k)+"}", coefficients[k]) for k in range(-4, 3)]
    return "\\begin{align*}\n" + "\\\\[2pt]\n".join(lines) + "\n\\end{align*}\n"


def obj(x):
    return {"sympy": str(x), "latex": tex_expr(x)}


def main():
    print("Reading bound sources", flush=True)
    current = read("benchmarks/self-energy/current.json")
    manifest = read("benchmarks/self-energy/"+current["manifest"])
    truth("current manifest SHA256", bindings["benchmarks/self-energy/"+current["manifest"]] == current["sha256"])
    targets = ("2e-Pi1-01", "2e-Pi1-02")
    records = {}
    for entry in manifest["entries"]:
        oid = entry["orbit_id"]
        if oid not in targets:
            continue
        for kind in ("raw", "prepared", "result"):
            record = entry[kind]
            data = read(record["path"])
            truth(oid+" "+kind+" SHA256", bindings[record["path"]] == record["sha256"])
            if kind == "result":
                records[oid] = data
    for rel in ("benchmarks/self-energy/conventions.md", "benchmarks/self-energy/report/common.tex",
                "benchmarks/self-energy/report/pi1-chain.tex", "benchmarks/renormalization/two-loop/derivation.md",
                "benchmarks/2ne/families/single-opposite/three-loop-r004.md"):
        bind(rel)
    existing_sum = read("benchmarks/renormalization/two-loop/result-r002.json")
    prepared=read("benchmarks/self-energy/pi1/prepared/2e-Pi1-r002.json")
    u,v=s.symbols("u v")
    fA=z*(1-u)*v
    fO=(1-u)*(1-v)+z*u*v
    for orbit,f,d,expected in zip(prepared["orbits"],(fA,fO),(z,-(1+z)),
            (2*z*fA,-2*(1+z)*fO+(1-2*e)*z)):
        oid=orbit["orbit_id"]
        numerator=2*(1+e)*d*f+f*s.diff(f,u,v)+(2*e-1)*s.diff(f,u)*s.diff(f,v)
        check(oid+" mixed derivative to direct kernel",numerator,expected)
        check(oid+" prepared rational prefactor",
              s.sympify(orbit["C"]["components"]["rational"],locals=LOCALS),-1/(16*e*(1-2*e)))
        truth(oid+" source order",orbit["J"]["iterated_order_outer_to_inner"]==(["s31","s41"] if oid==targets[0] else ["s21","s41"]),"convention")

    # Z**2 is analytic. Its fifth order is needed because I_A starts at eps^-3.
    logZ2 = sum(2*s.zeta(k)*e**k/k for k in range(2, 6))
    Z2 = exp_series(logZ2, 5)
    print("Checking Gamma expansion", flush=True)
    check("Gamma logarithm through eps5", s.series(2*s.loggamma(1-e)-2*s.EulerGamma*e,e,0,6).removeO(), logZ2)
    check("Gamma exponential recurrence", s.series(s.exp(logZ2),e,0,6).removeO(), Z2)
    # Original centered generating density, with nu=2eps.
    H = M01+e*M11+e**2*(M21/s.Integer(2)+M03/s.Integer(6))+e**3*(M31+M13)/6
    W,Q=s.symbols("W Q")
    density = s.series(s.exp(e*W)*s.sinh(e*Q)/e,e,0,4).removeO()
    moments = {W**k*Q**j: val for k,j,val in [(0,1,M01),(1,1,M11),(2,1,M21),(0,3,M03),(3,1,M31),(1,3,M13)]}
    for p in range(4):
        actual = sum(c*moments[W**k*Q**j] for (k,j),c in s.Poly(s.expand(density).coeff(e,p),W,Q).terms())
        check("centered moment density eps"+str(p),actual,s.expand(H).coeff(e,p))
    sech = 1-e**2*ell**2/2
    rational = s.series(1/(1-2*e),e,0,6).removeO()
    rational_A = s.series(1/((1-2*e)*(1+e)),e,0,6).removeO()
    A_poly = s.expand(-Z2*rational_A/(32*e**3))
    O_poly = s.expand(-Z2*(rational*H*sech/(16*e)+rational_A/(32*e**2)))
    geometry={targets[0]:{p:canonical(A_poly.coeff(e,p)) for p in range(-4,3)},
              targets[1]:{p:canonical(O_poly.coeff(e,p)) for p in range(-4,3)}}
    print("Checking current full-epsilon expressions", flush=True)

    # Bind the formula expanded here to the exact source expressions.
    gammaZ2=s.exp(-2*s.EulerGamma*e)*s.gamma(1-e)**2
    IA=-gammaZ2/(32*e**3*(1-2*e)*(1+e))
    H_exact=(-z-1)*z**(-e)*s.meijerg(((0,0,2*e),()),((0,2*e,2*e),()),z)/(s.gamma(1-2*e)*s.gamma(1+2*e)**2)
    IO=-gammaZ2/(16*e*(1-2*e))*(H_exact/s.cosh(e*s.log(z))+1/(2*e*(1+e)))
    for oid,exact in zip(targets,(IA,IO)):
        print(oid+" source and stored coefficients", flush=True)
        check(oid+" full epsilon source",exact,s.sympify(records[oid]["I"]["sympy"],locals=LOCALS))
        truth(oid+" physical sign and one Euler factor",records[oid]["normalization"]["assembly_sign"]==1 and records[oid]["normalization"]["Euler_count"]==1,"convention")
        for p in range(-6,1):
            old=s.sympify(records[oid]["laurent"]["coefficients"][str(p)]["expression"],locals=LOCALS)
            check(oid+" stored scalar eps"+str(p),geometry[oid].get(p,s.S.Zero),old)

    # Check all old two-loop scalar subtotal entries independently.
    for p in range(-4,1):
        check("stored Pi1 subtotal eps"+str(p),sum(geometry[oid][p] for oid in targets),
              s.sympify(existing_sum["group_coefficients"]["Pi1_without_B1"][str(p)]["sympy"],locals=LOCALS))

    weights={"A":(10-6*e)/(3-2*e),"f":-8*(1-e)/(3-2*e)}
    colors={}
    for oid in targets:
        print(oid+" color convolution", flush=True)
        colors[oid]={}
        for channel,weight in weights.items():
            w=s.series(weight,e,0,6).removeO().expand()
            cs={p:s.expand(sum(w.coeff(e,p-j)*geometry[oid][j] for j in range(-4,p+1))) for p in range(-4,3)}
            colors[oid][channel]=cs
            # Independent rational-denominator recurrence. This avoids asking a
            # generic series engine to repeatedly differentiate long M polynomials.
            n0,n1=(10,-6) if channel=="A" else (-8,8)
            direct={-5:s.S.Zero}
            for p in range(-4,3):
                direct[p]=s.expand((n0*geometry[oid][p]+n1*geometry[oid].get(p-1,0)+2*direct[p-1])/3)
            for p in range(-4,3):
                check(oid+" "+channel+" direct color convolution eps"+str(p),cs[p],direct[p])
    old_color={
      targets[0]:{"A":{-4:0,-3:-s.Rational(5,48),-2:-s.Rational(1,9),-1:-s.Rational(35,108)-5*s.pi**2/288,0:-s.Rational(89,162)-s.pi**2/54-5*s.zeta(3)/72},
                  "f":{-4:0,-3:s.Rational(1,12),-2:s.Rational(1,18),-1:s.Rational(11,54)+s.pi**2/72,0:s.Rational(49,162)+s.pi**2/108+s.zeta(3)/18}},
      targets[1]:{"A":{-4:0,-3:0,-2:-s.Rational(5,48),-1:-5*M01/24-s.Rational(1,9),0:-31*M01/72-5*M11/24-s.Rational(35,108)-5*s.pi**2/288},
                  "f":{-4:0,-3:0,-2:s.Rational(1,12),-1:M01/6+s.Rational(1,18),0:5*M01/18+M11/6+s.Rational(11,54)+s.pi**2/72}}}
    for oid in targets:
        for channel in weights:
            for p in range(-4,1):
                check(oid+" "+channel+" published two-loop color eps"+str(p),colors[oid][channel][p],old_color[oid][channel][p])
    truth("negative control: repeated historical minus changes the current leading pole",
          -geometry[targets[0]][-3] != -s.Rational(1,32),"negative")
    truth("negative control: freezing the gluon color weight changes the finite term",
          s.expand(weights["A"].subs(e,0)*geometry[targets[0]][0]-old_color[targets[0]]["A"][0]) != 0,"negative")

    preamble = r"""% Generated by workflow/checks/two-loop-report-20261015/build_self_energy.py.
The adjacent-edge representative \texttt{2e-Pi1-01} and opposite-edge
representative \texttt{2e-Pi1-02} contain the same one-loop transverse
self-energy operator. Its full colour factor is
\begin{equation}
 C_F\mathcal B_1(\eps)=C_F\left[a(\eps)C_A+b(\eps)T_Fn_f\right],\qquad
 a=\frac{10-6\eps}{3-2\eps},\quad b=-\frac{8(1-\eps)}{3-2\eps}.
\end{equation}
The two scalar integrals below are shared by both colour channels. They already
use the current physical sign, the orbit multiplicities $m_A=4$, $m_O=2$,
and one Euler compensation. The gluon--ghost and quark contributions are
contained in $\mathcal B_1$.

\subsection{The explicit \texorpdfstring{$C,J$}{C,J} connection}
Set $\ell=\log z$, $Z=e^{-\gamma_E\eps}\Gamma(1-\eps)$ and
$\mathcal R_1=-Z^2/[16\eps(1-2\eps)]$. With
$S_{14}=(\pi(-2v_{14}+i0))^{2\eps}$, the prepared single-representative
prefactor is $C=\lambda_b^2 C_F\mathcal B_1 S_{14}\mathcal R_1$,
where $\lambda_b=\alpha_s^b/\pi$.
This is the prepared bare normalization, before the $\overline{\mathrm{MS}}$
coupling conversion and the restoration of the physical $X_s,X_t$ scales
used in the main text.
The ordered parameters are $s=s_{31},t=s_{41}$ for $A$ and
$s=s_{21},t=s_{41}$ for $O$, with the inner $t$ integration performed first:
\begin{align}
 J_A&=\int_0^1\dd s\int_0^1\dd t\,
 \frac{z}{1+\eps}\,[z(1-s)t]^{2\eps-1}
 =\frac{z^{2\eps}}{4\eps^2(1+\eps)},\\
 J_O&=\int_0^1\dd s\int_0^1\dd t\,
 \frac{-2(1+z)f+(1-2\eps)z}{2(1+\eps)}f^{2\eps-2}\\
 &=z^\eps\mathcal H(2\eps,z)
   +\frac{1+z^{2\eps}}{4\eps(1+\eps)}.
\end{align}
Here $f=(1-s)(1-t)+zst$. In the last line the mixed derivative has been integrated to its four corner
values at $\Re\eps>0$, before meromorphic continuation. The symmetric mother
function and its centered moments are
\begin{align}
 \mathcal H(\nu,z)&=-(1+z)z^{-\nu/2}
  \int_0^1\dd s\int_0^1\dd t\,f^{\nu-1},\\
 \mathcal H(2\eps,z)&=M_{0,1}+\eps M_{1,1}
 +\eps^2\left(\frac{M_{2,1}}2+\frac{M_{0,3}}6\right)\notag\\
 &\quad+\frac{\eps^3}{6}(M_{3,1}+M_{1,3})+O(\eps^4).
\end{align}
These are the same $M$ functions defined in the common function dictionary.
The crossing combination gives
\begin{align}
 I_r&=\frac{m_r}{2}\mathcal R_1
 \frac{J_r(z)+z^{2\eps}J_r(1/z)}{1+z^{2\eps}},\\
 I_A&=-\frac{Z^2}{32\eps^3(1-2\eps)(1+\eps)},\\
 I_O&=-\frac{Z^2}{16\eps(1-2\eps)}
 \left[\frac{\mathcal H(2\eps,z)}{\cosh(\eps\ell)}
       +\frac1{2\eps(1+\eps)}\right].
\end{align}
Thus the orbit contribution is
$\lambda_b^2 C_F\mathcal B_1(S_{14}+S_{12})I_r$, with
$S_{12}=z^{2\eps}S_{14}$. No further rotation factor is required.

\clearpage
\subsection{Scalar results through \texorpdfstring{$\eps^2$}{epsilon squared}}
Write $I_r=\sum_{k=-4}^{2}i_{r,k}\eps^k+O(\eps^3)$.
The formulas use $M_{0,3}=-M_{0,1}^2$ and
$\ell^2=-2M_{0,1}-\pi^2$ to keep the coefficients in the original moment basis.
For the adjacent representative, $i_{A,k}$ is independent of $z$:
"""
    text=preamble+coeff_block("i_{A}",geometry[targets[0]])
    text+=r"For the opposite representative:"+"\n"+coeff_block("i_{O}",geometry[targets[1]])
    text+=r"""
\clearpage
\subsection{Colour-resolved results through \texorpdfstring{$\eps^2$}{epsilon squared}}
Define $aI_r=\sum c^{(A)}_{r,k}\eps^k$ and
$bI_r=\sum c^{(f)}_{r,k}\eps^k$, with both remainders $O(\eps^3)$.
These coefficients multiply $C_FC_A$ and $C_FT_Fn_f$, respectively.
The rational weights are expanded through $\eps^5$ for $A$ and through
$\eps^4$ for $O$ before truncation. Their positive powers multiply the
integral poles; setting $a,b$ to their four-dimensional values would alter
the answer.
"""
    for oid,short in zip(targets,("A","O")):
        if short=="O":
            text+="\\clearpage\n"
        text+="\\paragraph{"+("Adjacent" if short=="A" else "Opposite")+" representative \\texttt{"+oid+"}.}\n"
        for channel in weights:
            symbol="c^{("+channel+")}_{"+short+"}"
            text+=coeff_block(symbol,colors[oid][channel])
    # Combine the two subscripts generated by the generic renderer into one.
    for stem in ("i_{A}","i_{O}","c^{(A)}_{A}","c^{(f)}_{A}","c^{(A)}_{O}","c^{(f)}_{O}"):
        for p in range(-4,3):
            text=text.replace(stem+"_{"+str(p)+"}",stem[:-1]+","+str(p)+"}")
    TEX.write_text(text,encoding="utf-8")
    print("Writing self-energy artifacts", flush=True)
    data={"format":"wilson-two-loop-preread-self-energy/1","through_power":2,
          "remainder":"O(eps^3), fixed z>0","sign_convention":current["sign_convention"],
          "scalar_coefficients":{oid:{str(p):obj(x) for p,x in cs.items()} for oid,cs in geometry.items()},
          "color_coefficients":{oid:{ch:{str(p):obj(x) for p,x in cs.items()} for ch,cs in channels.items()} for oid,channels in colors.items()},
          "full_epsilon_formulae":{targets[0]:str(IA),targets[1]:str(IO)},
          "dependencies":bindings,"formal":"not_run"}
    (OUT/"self-energy-data.json").write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    bind("workflow/checks/two-loop-report-20261015/build_self_energy.py")
    report={"format":"wilson-report-algebra-checks/1","passed":all(c["passed"] for c in checks),
            "check_count":len(checks),"checks":checks,"dependencies":bindings,
            "limitations":["The full-epsilon integral source is reused, not independently recomputed.",
                           "The original M03 and log-squared identities are reused from the source dictionary.",
                           "No numerical integration or Lean proof is run by this builder."],"formal":"not_run"}
    (OUT/"self-energy-checks.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"passed":report["passed"],"checks":len(checks),"tex":str(TEX),"coefficient_count":42}))


if __name__=="__main__":
    main()
