"""Versioned self-energy sign normalization at the Wilson import boundary.

The upstream montopar exports are immutable evidence, not code owned by Wilson.
Current project inputs use i-included 1PI insertions with positive Dyson chaining.
The default 47 pure inputs and their preparation entry point are unchanged.
"""
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/self-energy'
CONVENTION='wilson-self-energy-1pi-insertion/2.0'
SPECS={
    '2e-Pi1':dict(folder='pi1',loops=2,count=2,insertions=1,operator='Pi1',kind='1pi',orders=[1],raw='2e-Pi1-r001.json'),
    '4e-Pi1':dict(folder='pi1',loops=3,count=11,insertions=1,operator='Pi1',kind='1pi',orders=[1],raw='4e-Pi1-r001.json'),
    '2e-Pi1+Pi1':dict(folder='pi1-chain',loops=3,count=2,insertions=2,operator='Pi1+Pi1',kind='one_particle_reducible',orders=[1,1],raw='montopar-2e-Pi1+Pi1-r001.json'),
    '2e-Pi2':dict(folder='pi2',loops=3,count=2,insertions=1,operator='Pi2',kind='1pi',orders=[2],raw='montopar-2e-Pi2-r001.json'),
}

def encoded(x):return (json.dumps(x,ensure_ascii=False,indent=2)+'\n').encode('utf8')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def require(ok,message):
    if not ok:raise ValueError(message)
def write_revision(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);blob=encoded(x)
    if p.exists():require(p.read_bytes()==blob,'immutable_revision_mismatch: '+str(p))
    else:p.write_bytes(blob)

def conversion_factor(family):
    require(family in SPECS,'unsupported_self_energy_family')
    return (-1)**SPECS[family]['insertions']

def prefactor_R(family):
    n=SPECS[family]['loops']
    if family.endswith('Pi1'):
        return f'(-1)^({n}-1)*G^{n}/(2^({n}+2)*eps*(1-2*eps))'
    sign='-' if family=='2e-Pi1+Pi1' else ''
    return sign+'G^3*U/(128*eps^2*(1-2*eps)^2)'

def _flip_C(C,factor):
    out=deepcopy(C)
    require(C.get('component_semantics')=='multiplicative','nonmultiplicative_C')
    order=C['factor_order'];parts=C['components']
    require(len(order)==len(set(order)) and set(order)==set(parts),'C_factor_order_mismatch')
    require(C['combined']=='*'.join('('+parts[k]+')' for k in order),'C_combined_mismatch')
    if factor==-1:out['components']['rational']='-('+parts['rational']+')'
    out['combined']='*'.join('('+out['components'][k]+')' for k in order)
    ledger=out['phase_ledger'];ledger['source_external_sign']='+1'
    for key in ['coordinate_product_before_extracting_vref','overall_dimensionless_sign',
                'momentum_prefactor_phase_relative_to_P2','coordinate_sign_before_scale_relative_to_P2',
                'overall_dimensionless_sign_relative_to_P2']:
        if factor==-1 and key in ledger:
            value=s.sympify(ledger[key].replace('^','**'))
            ledger[key]=str(-value)
    return out

def normalize_orbit(orbit,family,prepared):
    """Normalize a recognized legacy export exactly once; never guess its sign."""
    require(family in SPECS,'unsupported_self_energy_family');spec=SPECS[family]
    require('sign_convention' not in orbit,'self_energy_already_normalized')
    op=orbit['operator'];norm=orbit['normalization'];top=orbit['topology']
    require(top['emission_type']==family and top['nloops']==spec['loops'],'family_or_loop_mismatch')
    orders=op.get('subgraph_loop_orders',[op.get('subgraph_loop_order')])
    require(op['name']==spec['operator'] and op['kind']==spec['kind'] and orders==spec['orders']
            and op['insertion_count']==spec['insertions'],'operator_identity_mismatch')
    require(op['renormalization']=='bare' and op['gauge']=='Feynman' and op['dimension']=='4-2*eps','operator_scheme_mismatch')
    factor=conversion_factor(family)
    require(int(orbit['C']['phase_ledger']['source_external_sign'])==factor,'legacy_external_sign_mismatch')
    require(norm['euler_compensation_count']==int(prepared),'Euler_stage_mismatch')
    require(orbit['C']['components'].get('euler')==(f'exp(-{spec["loops"]}*EulerGamma*eps)' if prepared else None),'Euler_factor_mismatch')
    out=deepcopy(orbit);out['C']=_flip_C(orbit['C'],factor)
    if 'original_C' in out:out['original_C']=_flip_C(orbit['original_C'],factor)
    if 'original_normalization' in out:del out['original_normalization']
    out['sign_convention']=CONVENTION
    out['normalization']['scheme']=CONVENTION
    out['normalization']['assembly_sign']=1
    if 'shared_I_convention' in out['normalization']:
        out['normalization']['shared_I_convention']['common_R']=prefactor_R(family)
        out['normalization']['shared_I_convention']['status']='current_physical_I'
    if prepared:
        out['I_target']=dict(sym=orbit['c4_symmetry']['orbit_multiplicity'],
            a1=orbit.get('color_preparation',{}).get('a1','1'),total_loop_order=spec['loops'],
            R=prefactor_R(family),definition='(m*a1/2)*R*(J(z)+z^(n*eps)*J(1/z))/(1+z^(n*eps))',
            external_normalization='unchanged operator/color factors; see conventions.md',
            assembly_sign=1,status='derived_from_legacy_by_linearity',formal='not_run')
    out['convention_migration']=dict(from_external_sign=factor,to_external_sign=1,
        C_multiplier=factor,J_multiplier=1,I_multiplier=factor,
        operator_scalar_changed=False,color_changed=False,Euler_changed=False)
    return out

def _scaled_math(record,factor):
    out=deepcopy(record)
    if factor==1:return out
    require(any(k in out for k in ['expression','sympy','sympy_template']),'missing_math_expression')
    if 'expression' in out:
        out['expression']='-('+out['expression']+')'
    for key in ['sympy','sympy_template','classical_polylog_expression']:
        if key in out:out[key]='-('+out[key]+')'
    # Generate readable coefficients rather than double-negated printed expressions.
    if out.get('language')=='SymPy/str':out['expression']=str(s.expand(s.sympify(out['expression'])))
    tex_source=out.get('sympy',out.get('sympy_template'))
    if out.get('language')=='SymPy/str':tex_source=out['expression']
    if tex_source is not None:out['latex']=s.latex(s.sympify(tex_source))
    elif 'latex' in out:out['latex']='-\\left('+out['latex']+'\\right)'
    if 'hlog_terms' in out:
        for term in out['hlog_terms']:term['coefficient']=str(-Fraction(term['coefficient']))
    return out

def _scaled_laurent(block,factor):
    out=deepcopy(block)
    out['coefficients']={k:_scaled_math(v,factor) for k,v in block['coefficients'].items()}
    return out

def normalize_result(result,family,input_reference,legacy_reference):
    require('sign_convention' not in result,'result_already_normalized')
    require(result['revision']=='r001','unsupported_legacy_revision')
    factor=conversion_factor(family);out=deepcopy(result)
    for key in ['I','I_geom']:
        if key in out:out[key]=_scaled_math(out[key],factor)
    for key in ['laurent','I_geom_laurent']:
        if key in out:out[key]=_scaled_laurent(out[key],factor)
    for channel in out.get('color_channels',[]):
        if 'full_epsilon' in channel:channel['full_epsilon']=_scaled_math(channel['full_epsilon'],factor)
        if 'laurent' in channel:channel['laurent']=_scaled_laurent(channel['laurent'],factor)
        if 'coefficients' in channel:
            channel['coefficients']={k:_scaled_math(v,factor) for k,v in channel['coefficients'].items()}
    if 'verification' in out:out['historical_verification']=out.pop('verification')
    out['legacy_input']=out.pop('input')
    out['input']=input_reference;out['revision']='r002';out['sign_convention']=CONVENTION
    out['normalization']['scheme']=CONVENTION
    out['normalization']['assembly_sign']=1
    out['normalization']['R']=prefactor_R(family)
    out['normalization']['I']='(m*a1/2)*R*(J(z)+z^(n*eps)*J(1/z))/(1+z^(n*eps))'
    if family=='2e-Pi2':out['normalization']['shared_geometry']='I=I_geom shared by AA/Af/Ff; I_pi2=-I_chain; J kernels coincide'
    out['migration']=dict(legacy_result=legacy_reference,factor=factor,
        reason='linear C*J and unchanged external color/operator normalization',
        new_integration_performed=False,old_formal_status_not_inherited=True)
    out['verification']=dict(exact_migration='workflow/checks/self-energy-signs/checks.json',formal='not_run')
    out['formal']='not_run'
    return out

def publish(root=ROOT):
    root=Path(root);base=root/'benchmarks/self-energy';records=[];inputs={};published={}
    def reference(p):
        return dict(path=p.relative_to(root).as_posix(),sha256=digest(p))
    for family,spec in SPECS.items():
        folder=base/spec['folder'];rawold=folder/'sources'/spec['raw'];prepold=folder/'prepared'/(family+'.json')
        packages=[]
        for old,stage,prepared in [(rawold,'sources',False),(prepold,'prepared',True)]:
            src=read(old);inputs[old.relative_to(root).as_posix()]=digest(old)
            require(src['emission_type']==family and len(src['orbits'])==spec['count'],'source_package_mismatch')
            out=deepcopy(src);out['orbits']=[normalize_orbit(o,family,prepared) for o in src['orbits']]
            out['sign_convention']=CONVENTION;out['revision']='r002'
            out['normalization_contract']='../../conventions.md'
            out['legacy_package']=reference(old)
            path=folder/stage/(family+'-r002.json');write_revision(path,out)
            packages.append(path);published[path.relative_to(root).as_posix()]=digest(path)
        rawpath,preppath=packages
        for orbit in read(preppath)['orbits']:
            oid=orbit['orbit_id'];old=folder/'results'/oid/'r001.json';src=read(old)
            inputs[old.relative_to(root).as_posix()]=digest(old)
            require(src['orbit_id']==oid,'result_orbit_mismatch')
            require(src['input']['package_sha256']==digest(prepold),'legacy_result_input_mismatch')
            ref=dict(package=preppath.relative_to(root).as_posix(),package_sha256=digest(preppath),
                sign_convention=CONVENTION,raw_package=rawpath.relative_to(root).as_posix(),raw_package_sha256=digest(rawpath))
            out=normalize_result(src,family,ref,reference(old));path=old.with_name('r002.json');write_revision(path,out)
            published[path.relative_to(root).as_posix()]=digest(path)
            records.append(dict(orbit_id=oid,family=family,raw=reference(rawpath),prepared=reference(preppath),
                result=reference(path),legacy_result=reference(old),legacy_to_current=conversion_factor(family),assembly_sign=1))
    manifest=dict(format='wilson-self-energy-current/2',sign_convention=CONVENTION,
        operator='Pi_L=-i*(k^2*g-k*k)*P_L; all scalar/color/internal signs retained',
        dyson='D=D0+D0*Pi*D0+D0*Pi*D0*Pi*D0+...',
        external_source_sign=1,entries=records,legacy_inputs=inputs,files=published,
        producer=dict(path='preprocessing/self_energy.py',sha256=digest(Path(__file__))),
        raw_provenance='project-local normalization of frozen montopar output; external repository not modified',
        default_47_inputs_unchanged=True,formal='not_run')
    write_revision(base/'current-r002.json',manifest)
    # A small mutable selector is the only current-version pointer.
    (base/'current.json').write_bytes(encoded(dict(manifest='current-r002.json',sha256=digest(base/'current-r002.json'),sign_convention=CONVENTION)))
    return manifest

def load_current_result(orbit_id,root=ROOT):
    root=Path(root);base=root/'benchmarks/self-energy';pointer=read(base/'current.json')
    require(pointer['sign_convention']==CONVENTION,'unsupported_current_convention')
    p=base/pointer['manifest'];require(p.resolve().is_relative_to(base.resolve()),'manifest_outside_corpus')
    require(digest(p)==pointer['sha256'],'current_manifest_hash_mismatch');manifest=read(p)
    matches=[x for x in manifest['entries'] if x['orbit_id']==orbit_id]
    require(len(matches)==1,'unknown_or_duplicate_self_energy_orbit');entry=matches[0]
    for key in ['raw','prepared','result']:
        rec=entry[key];path=root/rec['path'];require(path.resolve().is_relative_to(base.resolve()),'artifact_outside_corpus')
        require(digest(path)==rec['sha256'],'current_artifact_hash_mismatch: '+key)
    out=read(root/entry['result']['path'])
    require(out['sign_convention']==CONVENTION and out['normalization']['assembly_sign']==1,'current_result_sign_mismatch')
    return out

if __name__=='__main__':
    manifest=publish();print(json.dumps({'convention':CONVENTION,'orbits':len(manifest['entries']),
        'flipped':sum(x['legacy_to_current']==-1 for x in manifest['entries']),'external_modified':False}))
