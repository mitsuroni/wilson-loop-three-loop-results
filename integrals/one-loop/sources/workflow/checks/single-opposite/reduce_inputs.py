"""Derive single-opposite weights from each actual ordered root, not an orbit lookup."""
from pathlib import Path
import sys,json,itertools
import sympy as sp

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"workflow/checks/empty-family"))
from setup_workspaces import sha,exact_json_hash
from empty_math import to_sympy,eps,z,ell,R2,G,G_exact,R2_exact,word_integral,exp_series
from preprocessing.ordered import integral_ast
DATA=ROOT/"benchmarks/2ne"
FAMILY=DATA/"families/single-opposite"
AT,BT,AS,BS,BAT,BAS=sp.symbols("AT BT AS BS BAT BAS")
VARS=(AT,BT,AS,BS,BAT,BAS)
HPARAMS={"H0":(0,0,0,0),"H1":(0,0,1,0),"H2":(0,0,2,0),
         "HW":(0,0,1,1),"HD":(1,0,1,0),"HX":(1,0,0,1)}
H={k:sp.Symbol(k) for k in HPARAMS}
HB={k:sp.Symbol(k+"bar") if k in ("H1","H2","HD") else H[k] for k in HPARAMS}

def load_inputs():
    mp=DATA/"prepared/manifest.json"
    man=json.loads(mp.read_text(encoding="utf-8"))
    assert man["format"]=="wilson-prepared-input-manifest/0.3"
    rows=[]
    for pkg in man["packages"]:
        path=DATA/"prepared"/pkg["path"]
        assert sha(path)==pkg["sha256"]
        obj=json.loads(path.read_text(encoding="utf-8"))
        for o in obj["orbits"]:
            if o["integral_structure"]["direct_opposite_edge_links"]["channel_partition"]==[1]:
                rows.append({"orbit_id":o["orbit_id"],"loop_order":obj["loop_order"],"orbit":o,
                    "input":{"manifest":"benchmarks/2ne/prepared/manifest.json","manifest_sha256":sha(mp),
                             "package":path.relative_to(ROOT).as_posix(),"package_sha256":sha(path),
                             "orbit_sha256":exact_json_hash(o),"normalization_version":man["recipe"]}})
    assert len(rows)==18
    return rows

def lower(word,a,b,ba):
    return {"":1,"A":a,"B":1-b,"AA":a*a/2,"BB":(1-b)**2/2,
            "AB":a*(1-b)-ba,"BA":ba}[word]
def upper(word,a,b,ba):
    return {"":1,"A":1-a,"B":b,"AA":(1-a)**2/2,"BB":b*b/2,
            "AB":R2-a+ba,"BA":a+b-a*b-R2-ba}[word]
def reflect(poly):
    return sp.expand(poly.subs({AT:BT,BT:AT,AS:BS,BS:AS,
                              BAT:AT+BT-AT*BT-R2-BAT,
                              BAS:AS+BS-AS*BS-R2-BAS},simultaneous=True))
def moment_key(powers):
    for key,p in HPARAMS.items():
        if powers in (p,p[2:]+p[:2]):return H[key]
        ref=(p[1],p[0],p[3],p[2])
        if powers in (ref,ref[2:]+ref[:2]):return HB[key]
    raise ValueError(("unknown moment",powers))
def to_moments(poly):
    assert not poly.has(BAT,BAS),"BA must first cancel under crossing"
    return sp.expand(sum(c*moment_key(p) for p,c in sp.Poly(poly,AT,BT,AS,BS).terms()))

def reduce(row):
    o,n=row["orbit"],row["loop_order"];j=o["J"]
    integral_ast(j)
    assert o["dispatch"]=="pure" and len(j["variables"])==2*n
    variables={s:sp.Symbol(s,positive=True) for s in j["variables"]}
    letters={};props=[];scale=sp.Integer(1);opp=None
    for fn,raw in j["propagator_polynomials"].items():
        f=to_sympy(raw,variables)
        pair=[(a,b) for a,b in itertools.combinations(variables,2)
              if sp.expand(f-((1-variables[a])*(1-variables[b])+z*variables[a]*variables[b]))==0]
        if pair:
            assert opp is None and len(pair)==1
            opp=pair[0]
            for v in opp:
                assert v not in letters;letters[v]="O"
            props.append({"id":fn,"type":"opposite","variables":list(opp),"residual":"0"})
        else:
            choices=[(a,b,q) for a,b in itertools.permutations(variables,2) for q in (sp.Integer(1),z)
                     if sp.expand(f-q*(1-variables[b])*variables[a])==0]
            assert len(choices)==1
            a,b,q=choices[0]
            assert a not in letters and b not in letters
            letters[a],letters[b]="A","B";scale*=q
            props.append({"id":fn,"type":"adjacent","A":a,"B":b,"q":str(q),"residual":"0"})
    assert opp is not None and set(letters)==set(variables)
    assert sp.expand(to_sympy(j["numerator"],variables)+(1+z)*scale)==0
    r=int(sp.degree(scale,z));assert scale==z**r
    child={};tops=[]
    for v,(lo,hi) in j["bounds"].items():
        assert lo=="0"
        if hi=="1":tops.append(v)
        else:
            assert hi in variables and hi not in child
            child[hi]=v
    seen=set();chains=[];weight=sp.Integer(1);marked=0
    for top in sorted(tops):
        vs=[];v=top
        while True:
            assert v not in seen;seen.add(v);vs.append(v)
            if v not in child:break
            v=child[v]
        vs=list(reversed(vs));word="".join(letters[v] for v in vs)
        assert word.count("O")<=1
        if "O" in word:
            i=word.index("O")
            a,b,ba=(AT,BT,BAT) if vs[i]==opp[0] else (AS,BS,BAS)
            pre,post=word[:i],word[i+1:]
            factor=lower(pre,a,b,ba)*upper(post,a,b,ba)
            marked+=1
        else:
            factor=sp.cancel(eps**len(word)*word_integral(word))
        weight*=factor
        chains.append({"variables_low_to_high":vs,"word":word,"normalized_weight":str(factor)})
    assert len(seen)==2*n and marked==2
    weight=sp.expand(weight);delta=sp.Rational(2*r-(n-1),2)
    # Physical crossing removes the antisymmetric BA part for this batch.
    phase=sp.exp(delta*eps*ell)
    paired=sp.expand((phase*weight+phase**-1*reflect(weight))/2)
    numerator=to_moments(paired)
    a1=sp.Rational(o["color_preparation"]["a1"]);sym=o["c4_symmetry"]["orbit_multiplicity"]
    assert str(a1)==o["color_preparation"]["coordinates"][0]
    ca,cf=sp.symbols("CA CF",positive=True)
    assert sp.expand(to_sympy(o["color_preparation"]["projected_color"],["CA","CF"])-a1*ca**(n-1)*cf)==0
    ef="Exp[-EulerGamma*eps]" if n==1 else f"Exp[-{n}*EulerGamma*eps]"
    expected=(f"({o['color_preparation']['projected_color']}) * ({ef}) * ((4*Pi*alphaS)^{n})"
              f" * (Gamma[1-eps]^{n}/(4*Pi^(2-eps))^{n}) * ((-2*v14+I0)^({n}*eps)) * ((-2)^(-{n}))")
    assert o["C"]["combined"]==expected
    kappa=sp.Rational(sym)*a1/(2*(-2)**n)
    I=kappa*G**n* numerator/(eps**(2*n-2)*sp.cosh(n*eps*ell/2))
    return {"chains":chains,"letters":letters,"propagators":props,"opposite":list(opp),"r":r,"delta":delta,
            "weight":weight,"paired_weight":paired,"moment_numerator":numerator,"I_compact":I,
            "a1":a1,"sym":sym,"kappa":kappa}

def meijer(p):
    a,b,c,d=p
    return -(1+z)*z**(-eps/2)*sp.meijerg([eps,-a*eps,-c*eps],[],[0,(b+1)*eps,(d+1)*eps],[],z)/(
        sp.gamma(1-eps)*sp.gamma(1+(a+b+1)*eps)*sp.gamma(1+(c+d+1)*eps))
def exact_I(a):
    subs={H[k]:meijer(p) for k,p in HPARAMS.items()}
    subs.update({HB[k]:meijer(p).subs(z,1/z) for k,p in HPARAMS.items() if HB[k]!=H[k]})
    return a["I_compact"].subs(subs).subs({G:G_exact,R2:R2_exact,ell:sp.log(z)})

if __name__=="__main__":
    for row in load_inputs():
        a=reduce(row)
        print(row["orbit_id"],[c["word"] for c in a["chains"]],"r=",a["r"],"W=",a["weight"],"I=",a["I_compact"],flush=True)
