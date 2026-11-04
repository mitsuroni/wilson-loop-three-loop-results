"""Check new P5 values against fixed references and labelled source bindings."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,ast,hashlib,json,math
import sympy as sp
PACK=Path(__file__).resolve().parents[1]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ptr(v,p):
    for t in p.strip('/').split('/') if p else []:
        t=t.replace('~1','/').replace('~0','~');v=v[int(t)] if isinstance(v,list) else v[t]
    return v
def vec(row):
    if row==0:return {}
    if isinstance(row,dict) and 'hlog_terms' not in row and row.get('language')=='SymPy/str':
        rho=sp.Symbol('rho');z=sp.Symbol('z',positive=True)
        def visit(n):
            if isinstance(n,ast.Constant) and type(n.value) is int:return sp.Integer(n.value)
            if isinstance(n,ast.Name) and n.id in ['pi','z']:return {'pi':sp.pi,'z':z}[n.id]
            if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.UAdd,ast.USub)):return visit(n.operand) if isinstance(n.op,ast.UAdd) else -visit(n.operand)
            if isinstance(n,ast.BinOp):
                a,b=visit(n.left),visit(n.right)
                if isinstance(n.op,ast.Add):return a+b
                if isinstance(n.op,ast.Sub):return a-b
                if isinstance(n.op,ast.Mult):return a*b
                if isinstance(n.op,ast.Div):return a/b
                if isinstance(n.op,ast.Pow):return a**b
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and len(n.args)==1 and not n.keywords:
                a=visit(n.args[0])
                if n.func.id=='zeta':return sp.zeta(a)
                if n.func.id=='log' and a==z:return rho
            raise ValueError('Unsupported scalar reference expression')
        expr=visit(ast.parse(row['expression'],mode='eval').body);out={}
        for powers,coefficient in sp.Poly(sp.expand(expr),rho,sp.pi,sp.zeta(3),sp.zeta(5)).terms():
            degree,pi_degree,z3,z5=powers
            if pi_degree%2:raise ValueError('Unexpected odd pi period')
            k=(0,)*degree,(pi_degree//2,z3,z5)
            out[k]=out.get(k,Q(0))+Q(str(coefficient))*math.factorial(degree)*6**(pi_degree//2)
        return {k:v for k,v in out.items() if v}
    if not isinstance(row,dict) or 'hlog_terms' not in row:raise ValueError('Expected complete explicit HLOG coefficient')
    out={}
    for t in row['hlog_terms']:
        k=tuple(t['word']),tuple(t['zeta_powers']);out[k]=out.get(k,Q(0))+Q(t['coefficient'])
        if any(a not in [0,-1] for a in k[0]):raise ValueError('Unexpected final alphabet')
    return {k:v for k,v in out.items() if v}
def key(c):
    if c.get('channel') in ['AA','Af','Ff','ff']:return c['channel']
    name=c.get('color_monomial',c.get('monomial',c.get('color',''))).replace(' ','').replace('**','^')
    aliases={'CA^2*CF':'AA','CA*CF*TF*nf':'Af','CF^2*TF*nf':'Ff','CF*(TF*nf)^2':'ff','CF*TF^2*nf^2':'ff'}
    if name not in aliases:raise ValueError('Unknown external colour '+name)
    return aliases[name]
def channels(r):
    values=r.get('color_channels',[])
    if isinstance(values,dict):return values
    return {key(c):c.get('coefficients',c.get('laurent',{}).get('coefficients')) for c in values}
def coefficients(r):
    if 'laurent' in r:return r['laurent']['coefficients']
    if 'I_geom_laurent' in r:return r['I_geom_laurent']['coefficients']
    raise ValueError('Missing scalar/geometric Laurent object')
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--results',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);a=parser.parse_args()
    if a.out.exists():raise ValueError('Output must be new')
    rows=read(a.results);inputs=read(PACK/'inputs.json');own={};i_count=c_count=0
    for bind in inputs['orbits']:
        oid=bind['orbit_id'];new=rows[oid];old=read(PACK/f'references/results/{oid}.json')
        n,o=coefficients(new),coefficients(old);positions=[]
        for p in range(-6,1):
            actual=vec(n[str(p)]);expected=vec(o[str(p)])
            if actual!=expected:raise ValueError('I coefficient mismatch '+oid+':'+str(p))
            positions.append({'power':p,'exact_reference_match':True,'HLOG_coordinates':len(actual)});i_count+=1
        nc,oc=channels(new),channels(old);channel_checks=[]
        for ch,cs in oc.items():
            if ch not in nc:raise ValueError('Missing channel '+ch)
            for p in range(-6,1):
                if vec(nc[ch][str(p)])!=vec(cs[str(p)]):raise ValueError('External operator/color mismatch '+oid+ch+str(p))
                c_count+=1
            channel_checks.append({'channel':ch,'positions':7,'exact_reference_match':True})
        values={}
        for name in ['native_prepared','native_raw','graph']:
            ref=bind[name];path=PACK/ref['path']
            if sha(path)!=ref['sha256']:raise ValueError('Changed original source')
            values[name]=ptr(read(path),ref['field'])
        native,raw=values['native_prepared'],values['native_raw']
        if native['original_C']!=raw['C'] or native['J']!=raw['J']:raise ValueError('Native preparation source bridge')
        if native['normalization']['euler_compensation_count']!=1 or native['normalization']['assembly_sign']!=1:raise ValueError('Euler/sign normalization')
        if native['preparation']['native_sign_conversion_count']!=0:raise ValueError('Unexpected sign patch')
        if native['orbit_id']!=oid or raw['orbit_id']!=oid:raise ValueError('Source identity')
        own[oid]={'I_positions':positions,'color_comparisons':channel_checks,'source_identity_CJ_Euler':True,
                  'quantity':'I_geom plus external contributions' if oid.startswith('2e-Pi2') else 'connected I plus external contributions'}
    report={'format':'wilson-p5-exact-reference-check/1','status':'passed','orbits':len(own),'I_positions':i_count,
            'colored_positions_compared':c_count,'results_sha256':sha(a.results),'by_orbit':own,
            'scope':'Fixed-result exact HLOG/source comparison after candidate freeze. Independent source/operator checks have their own report. No new original CAS integration or formal proof.',
            'formal':'not_run'}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode())
    print(json.dumps({k:report[k] for k in ['status','orbits','I_positions','colored_positions_compared']}))
if __name__=='__main__':main()
