"""Publish the verified physical channel, preserving every frozen source/result."""
from pathlib import Path
import hashlib
import json
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
inputs=read(HERE/'inputs.json');bridge=read(HERE/'physical-bridge.json');numeric=read(HERE/'numeric-checks.json')
assert len(inputs['parts'])==31 and bridge['passed'] and len(bridge['checks'])==31 and not bridge['RG_used_to_derive_sign']
for rel,h in {**inputs['source_sha256'],**bridge['source_sha256']}.items():assert sha(ROOT/rel)==h,rel
for job in ['physical','compact','source_controls','direct_export']:
    d=read(HERE/(job+'-run.json'));assert d['exit_code']==0 and not d['has_error']
    assert d['source_sha256']==sha(HERE/(job+'.mpl'))
physical=(HERE/'physical.log').read_text();assert 'PHYSICAL_CHANNEL_COMPLETED checks20 passed20' in physical
assert 'COMPACT_MN_SHAPE_AND_FINITE_EXACT_CHECKS_PASS=2' in (HERE/'compact.log').read_text()
def table(p):
    rows={}
    for line in p.read_text().splitlines():
        key,n,expr=line.split('\t');rows[key,int(n)]=expr
    return rows
ph=table(HERE/'physical-coefficients.txt');compact=table(HERE/'compact-coefficients.txt');old=table(HERE/'assembly-coefficients.txt')
checks={key[6:]:ex for (key,n),ex in ph.items() if key.startswith('check_')};assert len(checks)==20 and all(v=='0' for v in checks.values())
assert old['post',-5]=='1/27' and ph['omit_Gamma',-5]!='0'
assert ph['post',-6]==ph['post',-5]=='0'
lines=(HERE/'source-crossing.txt').read_text().splitlines();assert len(lines)==217 and all(x.split('\t')[-1]=='0' for x in lines)
direct=(HERE/'direct-export-checks.txt').read_text().splitlines();assert len(direct)==6 and all(x.split('\t')[-1]=='0' for x in direct)
assert numeric['all_pass'] and len(numeric['rows'])==63 and all(x['pass'] for x in numeric['rows'])
for rel,h in numeric['source_sha256'].items():assert sha(ROOT/rel)==h,rel
artifacts=['prepare.py','inputs.json','source_inputs.mpl','physical_bridge.py','physical-bridge.json','generate_physical.py','assemble.mpl','assemble.log','assembly-coefficients.txt','physical.mpl','physical.log','physical-run.json','physical-coefficients.txt','compact.mpl','compact.log','compact-run.json','compact-coefficients.txt','source_controls.mpl','source_controls.log','source-crossing.txt','direct_export.mpl','direct_export.log','direct-export-checks.txt','numeric.py','numeric-checks.json','publish.py']
result={'format':'wilson-three-loop-physical-color-result/1','revision':'r001','date':'2026-10-06','channel':'CF*CA*TF*nf','status':'computational_complete','coupling':'lambda=alpha_s_MSbar/pi','domain':'z>0, rho=log(z), L=(Ls+Lt)/2; fixed source meromorphic convention','source_graph_count':31,'source_group_counts':{'Pi1':11,'chain':2,'Pi2':2,'legs':12,'Gamma':4},'normalization':'external lambda^3*CF*CA*TF*nf extracted; each I already includes C4/web and one Euler compensation; S3 restored and lower S2/S1 preserved','physical_source_weights':bridge['weights'],'source_order_correction':'TGV-leg and Gamma1 stored I use chronological three-slot fTrace; physical P trace reverses the three slots. Use multiplier -1 without modifying I. Derived independently of RG.','all_original_inputs_unchanged':True,'upstream_source_colour_order_repaired':False,'raw_unit_coefficients':{str(n):ph['raw',n] for n in range(-6,1)},'bare_full_scale_coefficients':{str(n):ph['bare',n] for n in range(-6,1)},'post_counterterm_coefficients':{str(n):ph['post',n] for n in range(-6,1)},'log_Z4_coefficients':{str(n):ph['post',n] for n in range(-6,0)},'finite_GPL':ph['finite',0],'shape_GPL':ph['shape',0],'finite_MN':compact['finite',0],'shape_MN':compact['shape',0],'M_N_semantics':'original fixed moments at the same z; N12=M13-rho*M12; M01=-(rho^2+pi^2)/2','k3':ph['k3',0],'h3':ph['h3',0],'gamma4_convention':'gamma4_3=2*h3','epsilon_remainder':'O(eps)','checks':checks,'verification':{'physical_bridge_exact_controls':31,'individual_crossing_positions':217,'complete_physical_exact_checks':20,'compact_back_substitutions':2,'direct_export_checks':6,'independent_numeric_comparisons':63,'negative_controls':['As-exported weights fail with post epsilon^-5=1/27','Omitting Gamma1 leaves a nonzero epsilon^-5']},'cusp_reference':{'url':'https://arxiv.org/html/1510.07803','equation':'134','scheme':'MS, alpha_s/pi','coefficient':'-209/216+5*pi^2/54-7*zeta(3)/6'},'numeric_samples_at_L0':numeric['samples_at_L0'],'source_sha256':inputs['source_sha256'],'physical_bridge_source_sha256':bridge['source_sha256'],'evidence_sha256':{n:sha(HERE/n) for n in artifacts},'derivation':'workspaces/CF-CA-TF-nf/derivation.md','evidence':'workflow/checks/w3-ca-nf/README.md','formal':'not_run','MCP':'unavailable in this session','limitations':['Current emitted graph set and adopted scalar inputs; no independent completeness theorem for the producer','General complex kinematics and strict global source regulator-removal bridge remain open','Numeric series refinement is not an interval certificate or a re-integration of all original diagrams','Upstream production colour order has not been rewritten; future physical assembly must retain this explicit source-to-physical bridge'],'full_w3_complete':False}
family=ROOT/'benchmarks/renormalization/three-loop';dest=family/'results/CF-CA-TF-nf/r001.json';dest.parent.mkdir(parents=True,exist_ok=True)
payload=(json.dumps(result,ensure_ascii=False,indent=2)+'\n').encode()
if dest.exists():assert dest.read_bytes()==payload,'immutable r001 already differs'
else:dest.write_bytes(payload)
current={'format':'wilson-three-loop-physical-colors-current/1','date':'2026-10-06','channels':{'CF*CA*TF*nf':{'result':str(dest.relative_to(family)).replace('\\','/'),'sha256':sha(dest),'status':'computational_complete'},'CF^2*TF*nf':{'result':'matter-channels-r001.json','field':'cff','sha256':sha(family/'matter-channels-r001.json'),'status':'previous_computational_completion'},'CF*(TF*nf)^2':{'result':'matter-channels-r001.json','field':'nf2','sha256':sha(family/'matter-channels-r001.json'),'status':'previous_computational_completion'},'CF*CA^2':{'status':'incomplete','missing_complete_integrals':10,'missing_family':'double-TGV'}},'complete_color_channels':3,'total_color_channels':4,'full_w3_complete':False,'formal':'not_run','physical_bridge':'workflow/checks/w3-ca-nf/physical-bridge.json'}
(family/'current.json').write_text(json.dumps(current,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'published':str(dest.relative_to(ROOT)),'sha256':sha(dest),'complete_color_channels':3,'cusp':result['k3'],'h3':result['h3'],'formal':'not_run'}))
