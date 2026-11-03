"""Check formula-preserving manuscript edits; do not rerun historic CAS/numerics."""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    data=read(HERE/'input.json')
    assert sha(ROOT/data['archive'])==data['archived_sha256']
    ref=ROOT/data['style_reference']
    if ref.exists():assert sha(ref)==data['style_reference_sha256']
    e,l,Z,Q,Jv,VM,VN,VM2,chi,Nraw=s.symbols('eps ell Z Q J VM VN VM2 chi Nraw')
    # Current C/I and the self-mirror scaling.
    initial=-Z**2/s.Integer(8)*(Jv+s.exp(2*e*l)*s.exp(-2*e*l)*Jv)/(1+s.exp(2*e*l))
    final=-Z**2*Q/(8*s.cosh(e*l))
    assert s.simplify((initial.subs(Jv,s.exp(e*l)*Q)-final).rewrite(s.exp))==0
    assert s.Rational(2)*s.Rational(-1,2)/(2*(-2)**2)==s.Rational(-1,8)
    # V_N in the new exposition includes epsilon^-2; old script Nraw did not.
    oldQ=VM**2/2-2*Nraw/e**2+2*chi*VM2/e
    newQ=VM**2/2-2*VN+2*chi*VM2/e
    assert s.expand(oldQ.subs(Nraw,e**2*VN)-newQ)==0
    assert s.expand(oldQ.subs(Nraw,e**2*VN)-(VM**2/2-2*VN+chi*VM2/e))!=0
    # One-loop hyperbolic decomposition.
    logb,logbc=s.symbols('logbeta logbetabar',real=True)
    numerator=s.exp(-e*l/2)*(s.exp(e*(l+logbc))-s.exp(e*logb))
    hyper=2*s.exp(e*(logb+logbc)/2)*s.sinh(e*(l+logbc-logb)/2)
    assert s.simplify(s.expand((numerator-hyper).rewrite(s.exp)))==0
    # Central-variable geometry and the explicitly integrable rearrangement.
    p,q,z=s.symbols('p q z')
    u,v=(p+1)/(1+z),(q+1)/(1+z)
    assert s.factor((1-u)*(1-v)+z*u*v-(z+p*q)/(1+z))==0
    aa,bb,bq,cc=s.symbols('A B Bq C')
    original=(cc-bq)*(aa-cc)
    rearranged=(aa+bb)*(cc-1)-(cc**2-1)-(aa-1)*(bq-1)
    assert s.expand(rearranged-original-(bb-bq)*(cc-1))==0
    # Boundary ledger of the 4D -> 2D -> ordered-Beta evaluation.
    jm,j0,jp,Fsum,Tsum,Tone,Tdouble=s.symbols('jm j0 jp Fsum Tsum Tone Tdouble')
    Ksum=2*Fsum-2*Tsum-j0*(jm+jp)
    ledger=Ksum+(jm+jp)**2/2-j0*(jm+jp)+2*(Tsum-2*Tone)-Tdouble
    PsiDiag=Fsum-2*(Tone+j0**2/2)
    compact=(jm+jp-2*j0)**2/2+2*PsiDiag-Tdouble
    assert s.expand(ledger-compact)==0
    # Examples are generated with the new regular V_M/V_N definitions.
    M,N=s.Function('M'),s.Function('N')
    z2,z3=s.zeta(2),s.zeta(3)
    vseries=M(0,1)+e*M(1,1)/2
    nseries=N(0,2)/2+e*N(1,2)/2
    qseries=s.expand(vseries**2/2-2*nseries+2*(2*z2-4*z3*e)*vseries.subs(e,2*e))
    c0=s.expand((-qseries/8).coeff(e,0).subs(N(0,2),M(0,1)**2/3+4*z2*M(0,1)))
    c1=s.expand((-qseries/8).coeff(e,1))
    assert s.expand(c0+M(0,1)**2/48)==0
    target1=-M(0,1)*M(1,1)/16+z3*M(0,1)-s.pi**2*M(1,1)/12+N(1,2)/8
    assert s.simplify(c1-target1)==0
    prior=read(ROOT/'workflow/checks/4e04-direct-mn/validation.json')
    assert prior['passed']
    text=(ROOT/data['target']).read_text(encoding='utf8')
    required=['原四重有序积分','4D → 2D','对排序积分求z导数','N自然出现',
              '端点项χ','完整母函数','通项公式','附录A','附录B','4e04-direct-mn-r001.md']
    assert all(token in text for token in required)
    assert '\\tag{R' not in text
    report={'passed':True,'rewrite_only':True,'new_CAS_or_numerical_integrations':False,
      'archive_sha256':sha(ROOT/data['archive']),'style_reference_changed':False,
      'normalization_equivalent':True,'regular_N_generator_scaling_equivalent':True,
      'one_loop_hyperbolic_identity':True,'twoD_regularization_and_Beta_boundary_ledger':True,
      'first_two_I_coefficients':{'C0':str(c0),'C1':str(c1)},
      'dropped_boundary_factor_negative_control':True,
      'prior_math_validation_reused':'workflow/checks/4e04-direct-mn/validation.json',
      'formal':'not_run'}
    (HERE/'checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('presentation identities and archived evidence passed; no new integral computation')

if __name__=='__main__':main()
