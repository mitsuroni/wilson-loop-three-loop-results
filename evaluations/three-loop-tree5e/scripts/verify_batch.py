"""Compare fixed P6 candidate HLOG coefficients and all source bindings exactly."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,math
import sympy as sp
PACK=Path(__file__).resolve().parents[1]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ptr(v,p):
    for t in p.strip('/').split('/') if p else []:
        t=t.replace('~1','/').replace('~0','~');v=v[int(t)] if isinstance(v,list) else v[t]
    return v
def vector(terms):
    if isinstance(terms,dict):terms=terms['hlog_terms']
    out={}
    for t in terms:
        word=tuple(Q(str(x)) for x in t.get('word',t.get('Hlog_word',[])))
        if 'period_powers' in t:periods=tuple(t['period_powers']);factor=Q(1)
        elif 'zeta_powers' in t:periods=tuple(t['zeta_powers'])+(0,0,0);factor=Q(1)
        else:
            degree=t.get('pi_power',0)
            assert degree%2==0
            periods=(degree//2,t.get('zeta3_power',0),t.get('zeta5_power',0),0,0,0);factor=Q(6)**(degree//2)
        assert len(periods)==6
        key=word,periods;out[key]=out.get(key,Q(0))+Q(t['coefficient'])*factor
    return {k:v for k,v in out.items() if v}
def reference_coefficients(r):
    if 'hlog_coefficients' in r:return r['hlog_coefficients']
    return r['laurent']['coefficients']
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise ValueError('Use a new output file')
    rows=read(a.results);own={};count=0;alphabet=set();extra=[]
    for b in read(PACK/'inputs.json')['orbits']:
        oid=b['orbit_id'];old=read(PACK/f'references/results/{oid}.json');new=rows[oid];prior=reference_coefficients(old);cs=new['laurent']['coefficients'];positions=[]
        for power in range(-6,1):
            x,y=vector(cs[str(power)]),vector(prior[str(power)])
            if x!=y:
                keys=set(x)|set(y);res=[(str(k),str(x.get(k,0)-y.get(k,0))) for k in keys if x.get(k,0)!=y.get(k,0)]
                raise ValueError('Coefficient residual '+oid+':'+str(power)+' '+str(res[:5]))
            for w,_ in x:alphabet.update(w)
            if any(set(w)-{0,-1} or any(ps[3:]) for w,ps in x):extra.append({'orbit_id':oid,'power':power})
            positions.append({'power':power,'exact_reference_match':True,'HLOG_coordinates':len(x)});count+=1
        vals={}
        for name in ('native_prepared','native_raw','graph'):
            ref=b[name];path=PACK/ref['path'];assert sha(path)==ref['sha256'];vals[name]=ptr(read(path),ref['field'])
        native,raw=vals['native_prepared'],vals['native_raw']
        assert native['orbit_id']==raw['orbit_id']==oid
        assert native['original_C']==raw['C'] and native['J']==raw['J']
        assert native['normalization']['assembly_sign']==1 and native['normalization']['euler_compensation_count']==1
        assert native['preparation']['native_sign_conversion_count']==0
        own[oid]={'I_positions':positions,'source_identity_CJ_Euler':True,'quantity':'physical connected I; not ordered J or an IBP block'}
    assert len(own)==44 and count==308
    report={'status':'passed','orbits':44,'I_positions':count,'by_orbit':own,'alphabet':[str(x) for x in sorted(alphabet)],
        'extra_letter_or_period_positions':extra,'result_sha256':sha(a.results),
        'scope':'Fixed-current reference comparison after candidate freeze and byte/source identity. Representation, source controls and old numerical provenance are separate.', 'formal':'not_run'}
    a.out.write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode());print(json.dumps({k:report[k] for k in ('status','orbits','I_positions','extra_letter_or_period_positions')}))
if __name__=='__main__':main()
