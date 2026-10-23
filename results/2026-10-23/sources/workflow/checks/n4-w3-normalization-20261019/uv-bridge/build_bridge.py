"""Derive the coordinate-scale bridge from bare factors before comparing any h coefficient."""
from pathlib import Path
import json,hashlib,gzip,importlib.util,re
import sympy as s
import mpmath as mp
H=Path(__file__).resolve().parent;R=H.parents[3];P=R/'.wilson/literature/n4-uv-bridge-20261019'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
bindings={}
def rd(p):
 bindings[p.relative_to(R).as_posix()]=sha(p);return json.loads(p.read_text(encoding='utf-8'))
spec=importlib.util.spec_from_file_location('matter',R/'workflow/checks/matter-rg-20261017/check.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.ROOT=R
bindings['workflow/checks/matter-rg-20261017/check.py']=sha(R/'workflow/checks/matter-rg-20261017/check.py')
e,L,rho,CA,CF,T=m.eps,m.L,m.rho,m.CA,m.CF,m.T
z2,z3,z4,z5=s.zeta(2),s.zeta(3),s.zeta(4),s.zeta(5)
one=rd(R/'benchmarks/renormalization/one-loop/result-r001.json')
two=rd(R/'benchmarks/renormalization/two-loop/result-r002.json')
n4two=rd(R/'benchmarks/renormalization/n4-two-loop/comparison-r002.json')
delta=rd(R/'workflow/checks/n4-five-family-20261019/renormalized-difference-r001.json')
checks=[]
def eq(name,value):
 value=m.norm(s.expand(value));assert value==0,(name,str(value))
 checks.append({'name':name,'passed':True,'residual':str(value)})
# Real source bytes and exact text anchors, including the old single-TeX gzip format.
original=gzip.decompress((P/'9303314-source.tar').read_bytes())
assert original==(P/'9303314-source.tex').read_bytes()
tex=original.decode('latin1')
assert r'\def\MS{\overline{MS}}' in tex and r'\frac{g^2}{(2\pi)^{D/2}}' in tex and r'\lab{res:1}' in tex
fgm=(P/'1909-asset-36.tex').read_text(encoding='utf-8')
relations=(P/'1909-asset-37.tex').read_text(encoding='utf-8')
assert r'(4\pi e^{\gamma_E})^\epsilon' in fgm and r'\label{eq:w1bare}' in fgm
assert r'\label{eq:WBoxren}' in relations and r'\Gamma_\Box=4' in relations
sources={p.as_posix():sha(p) for p in [P/'9303314-source.tar',P/'9303314-source.tex',P/'1909-source.tar',P/'1909-asset-36.tex',P/'1909-asset-37.tex']}
# Coupling convention at the physical MSbar scale; all bases here are positive constants.
gE=s.EulerGamma
bare_coupling_log=gE-s.log(4*s.pi)
raw_logs={'project':s.log(s.pi),'closed':s.log(2*s.pi),'FGM_local':s.log(4*s.pi)}
scale_logs={'project':2*gE-s.log(4),'closed':2*gE-s.log(2),'FGM_local':2*gE}
for name in raw_logs:
 eq(name+' bare factor to native Euler Gamma',s.expand_log(raw_logs[name]+bare_coupling_log+gE-scale_logs[name],force=True))
eq('closed scale is twice project scale',scale_logs['closed']-scale_logs['project']-s.log(2))
eq('FGM local scale is twice original closed scale',scale_logs['FGM_local']-scale_logs['closed']-s.log(2))
# q=n.nstar, s=2q. Compare log arguments after the independently derived scale map.
q,M=s.symbols('q M',positive=True)
mux2=s.exp(2*gE)*M/4;muClosed2=s.exp(2*gE)*M/2;muF2=s.exp(2*gE)*M
eq('same invariant log, closed vs project',muClosed2*q-mux2*(2*q))
eq('same invariant log, FGM local vs project',muF2*(4*q)/4-2*muClosed2*q) # Q=4uv vs q bare units; scalar distance convention identity
# Primitive finite parallelogram at z=-1: opposite tangent contractions vanish.
ell=s.Symbol('ell')
G=s.exp(-gE*e)*s.gamma(1-e)
rawW=-G*(s.exp(e*ell)+s.exp(e*(ell+s.I*s.pi)))/e**2
series=s.expand(s.series(rawW,e,0,1).removeO())
F1=m.norm(m.parse(one['laurent_unit_coefficients']['0']).subs({s.Symbol('L_s'):L+rho/2,s.Symbol('L_t'):L-rho/2}))
eq('closed one-loop double pole',series.coeff(e,-2)+2)
eq('closed one-loop single pole',series.coeff(e,-1)+2*ell+s.I*s.pi)
eq('raw finite maps to accepted one loop at rho=i*pi',series.coeff(e,0).subs(ell,L-s.I*s.pi/2)-F1.subs(rho,s.I*s.pi))
# Two-loop accepted N4 values provide holdouts, not parameters to solve the scale map.
eq('accepted N4 two-loop k2',m.parse(n4two['cusp_lambda2_over_CFCA'])+s.pi**2/12)
eq('accepted N4 two-loop h2',m.parse(n4two['gamma4_lambda2_over_CFCA'])/2+7*z3/2)
# Published finite-line geometry relation Gamma_Box=4 f_eik and mu d/dmu = 2 d/dL + beta_mu d/dlambda.
a=s.Symbol('a')
B=12*z3*a**2-(80*z5+16*z2*z3)*a**3
G0=-4*z3*a**2+(32*z5+s.Rational(80,3)*z2*z3)*a**3
feik=s.expand(G0-2*B)
hN={2:m.norm(2*feik.coeff(a,2)/4**2),3:m.norm(2*feik.coeff(a,3)/4**3)}
eq('N4 h2 from independent G0 and B holdout',hN[2]+7*z3/2)
h3N=hN[3];k3N=s.Rational(11,8)*z4
# Explicit published QCD Gamma_Box/2 coefficients (1909.00697 eq5.15 and5.16).
h2ref=CA*(s.Rational(101,27)-s.Rational(11,12)*z2-7*z3/2)+T*(-s.Rational(28,27)+z2/3)
eq('two-loop both QCD colors physical geometry holdout',m.parse(two['h2'])/CF-h2ref)
href={
'AA':s.Rational(11,5)*z2**2+s.Rational(11,6)*z2*z3-s.Rational(6325,1296)*z2-s.Rational(329,24)*z3+6*z5+s.Rational(136781,23328),
'Af':-s.Rational(6,5)*z2**2+s.Rational(707,324)*z2+s.Rational(91,54)*z3-s.Rational(5921,5832),
'Ff':s.Rational(2,5)*z2**2+z2/4+s.Rational(19,9)*z3-s.Rational(1711,432),
'ff':-s.Rational(5,27)*z2+s.Rational(14,27)*z3-s.Rational(260,729)}
kref={
'AA':s.Rational(245,96)-s.Rational(67,36)*z2+s.Rational(11,24)*z3+s.Rational(11,8)*z4,
'Af':-s.Rational(209,216)+s.Rational(5,9)*z2-s.Rational(7,6)*z3,
'Ff':-s.Rational(55,48)+z3,'ff':-s.Rational(1,27)}
C3=s.Symbol('C3_N4')
Fn=-k3N*(L**2-rho**2/4)-h3N*L+C3
S3=s.exp(3*e*(L+rho/2))+s.exp(3*e*(L-rho/2))
Wn=S3*(-k3N/(9*e**2)-h3N/(6*e))+k3N*rho**2/2+C3
wnser=s.expand(s.series(Wn,e,0,1).removeO())
eq('N4 finite strict Laurent after full scale powers',wnser.coeff(e,0)-Fn)
total={}
for c,row in delta['channels'].items():
 kh=m.parse(row['delta_k3'])+(k3N if c=='AA' else 0)
 hh=m.parse(row['delta_h3'])+(h3N if c=='AA' else 0)
 eq(c+' independent full QCD cusp',kh-kref[c])
 eq(c+' independent full QCD Gamma_Box/2',hh-href[c])
 post={p:m.norm(m.parse(v)+(wnser.coeff(e,int(p)) if c=='AA' else 0)) for p,v in row['post_counterterm'].items()}
 total[c]={'post_counterterm':{p:str(v) for p,v in post.items()},'k3':str(m.norm(kh)),'h3':str(m.norm(hh))}
 # Independent scale dependence from full QCD lower-order result and published k3/h3.
 b0=(11*CA-4*T)/12;b1=(34*CA**2-20*CA*T-12*CF*T)/48
 extract={'AA':lambda x:s.expand(x).coeff(CA,2),'Af':lambda x:s.expand(x).coeff(CA).coeff(T),
          'Ff':lambda x:s.expand(x).coeff(CF).coeff(T),'ff':lambda x:s.expand(x).coeff(T,2)}[c]
 mpath=R/'workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt'
 bindings[mpath.relative_to(R).as_posix()]=sha(mpath)
 moment=m.norm(m.parse(mpath.read_text().split('=',1)[1].strip().rstrip(';')))
 F2=m.norm(m.parse(two['F2_over_CF']).subs({s.Symbol('M_1_1'):moment}));Fone=m.parse(two['F1_over_CF'])
 # Moment form is reused with registered exact M/N maps by delta, so compare derivative in same representation.
 res=s.diff(post['0'],L)-extract(s.expand(2*b0*F2+b1*Fone))+2*kref[c]*L+href[c]
 # Only elementary M01 can survive L differentiation.
 res=res.subs({s.Symbol('M01'):-(rho**2+s.pi**2)/2})
 eq(c+' full finite RG with independent references',res)
eq('C3 affects only AA finite constant',s.diff(s.sympify(total['AA']['post_counterterm']['0']),C3)-1)
# Numerical primitive cusp controls, independent of symbolic Laurent extraction.
mp.mp.dps=80;num=[]
for ep in map(mp.mpf,['0.2','0.35','0.55']):
 for qv in map(mp.mpf,['0.7','1.3']):
  Ms=mp.mpf('1.1')
  line=mp.quad(lambda t:mp.exp(-ep*t),[0,mp.inf])
  old=-mp.gamma(1-ep)*mp.exp(ep*(mp.euler-mp.log(2)))*(Ms*qv)**ep*(1+mp.exp(mp.j*mp.pi*ep))*line**2
  root=-mp.exp(-mp.euler*ep)*mp.gamma(1-ep)*(mp.exp(2*mp.euler)*Ms*qv/2)**ep*(1+mp.exp(mp.j*mp.pi*ep))/ep**2
  err=abs(old-root)/max(1,abs(root));assert err<mp.mpf('1e-65')
  num.append({'eps':str(ep),'q':str(qv),'normalized_error':mp.nstr(err,12),'passed':True})
result={'format':'wilson-N4-three-loop-UV-coordinate-bridge-research/1','revision':'r002','date':'2026-10-19','all_passed':True,
 'exact_checks':checks,'primitive_cusp_numeric_controls':num,'source_bindings':bindings,'downloaded_source_sha256':sources,
 'bare_coupling':'g_B^2 mu0^(2eps)=4 pi^2 lambda (mu_MSbar^2 exp(EulerGamma)/(4pi))^eps Z_lambda',
 'Euler_G':'exp(-EulerGamma eps) Gamma(1-eps), exactly one per loop',
 'coordinate_scales':{'mu_x2':'exp(2 gamma_E) mu_MSbar2/4','mu_closed2':'exp(2 gamma_E) mu_MSbar2/2','mu_FGM_local2':'exp(2 gamma_E) mu_MSbar2'},
 'log_connection':'s=2 n.nstar; mu_closed2*(n.nstar)=mu_x2*s; mu_closed2=mu_FGM_local2/2. Distinguish local integration scale from inherited closed-contour geometric scale.',
 'scale_map_not_fitted':'Obtained from raw pi, 2pi and 4pi factors and the common bare coupling before comparing any h coefficient. Full primitive one-loop finite branch and low-loop data are holdouts.',
 'geometry_source':{'url':'https://arxiv.org/pdf/1909.00697','equations':['(5.5)','(5.8)-(5.11)','(5.15)','(5.16)'],
 'operator':'closed four-finite-segment Wilson loop; Gamma_Box=4*(Gamma_cap-Gamma_wedge)=4*f_eik, not four independent infinite wedge constants'},
 'N4_h3_over_CFCA2':str(h3N),'N4_k3_over_CFCA2':str(m.norm(k3N)),
 'N4_F3_over_CFCA2':str(m.norm(Fn)),'h3_adopted_research':True,'C3_determined':False,
 'adoption_boundary':'Published finite-line geometry factorisation/UV-locality relation and independently known N4 G0 and B through three loops; checked against raw low-loop normalization and four QCD h3 channels. Not a primitive three-loop UV integration or formal proof.',
 'printed_1712_diagnostic':'Phase-1 literal-sign residual remains historical. This bridge uses the explicit finite-line RG object, not a silent sign flip in (2.7). No author erratum is asserted.',
 'formal':'not_run','selector_changed':False}
def store_json(name,value):
 p=H/name
 if p.exists():
  assert json.loads(p.read_text(encoding='utf-8'))==value,('immutable result differs',name)
 else:
  p.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
store_json('uv-coordinate-bridge-r002.json',result)
out={'format':'wilson-QCD-three-loop-total-up-to-N4-basepoint-research/1','revision':'r001','date':'2026-10-19',
 'status':'complete_poles_and_finite_shape_up_to_one_undetermined_constant_research','unknown':'C3_N4 in AA finite term only',
 'colors':delta['colors'],'channels':total,'bridge_sha256':sha(H/'uv-coordinate-bridge-r002.json'),'formal':'not_run',
 'absolute_F3_complete':False,'selector_changed':False}
store_json('qcd-total-up-to-C3-r001.json',out)
print(json.dumps({'all_passed':True,'exact_checks':len(checks),'numeric_controls':len(num),'h3_N4':str(h3N),'full_QCD_h3_channels':4,'C3_determined':False}))
