"""Read published ancillary files as data, fix conventions, compare analytically."""
from pathlib import Path
import hashlib
import json
import re
import requests
import sympy as s
from sympy.parsing.mathematica import parse_mathematica
from calculate import ROOT,HERE,DEST,data,encode,sha,obj

CACHE=ROOT/'.tmp/rg-two-loop-papers'
DOI='https://doi.org/10.5281/zenodo.17200473'

def fetch():
    CACHE.mkdir(parents=True,exist_ok=True)
    metapath=CACHE/'zenodo-17200473.json'
    if not metapath.exists():
        response=requests.get('https://zenodo.org/api/records/17200473',timeout=25)
        response.raise_for_status();metapath.write_bytes(response.content)
    meta=data(metapath);records=[]
    for item in meta['files']:
        if item['key'] not in ['WRen.m','ADs.m','README.txt']:continue
        path=CACHE/item['key']
        if not path.exists():
            response=requests.get(item['links']['self'],timeout=25);response.raise_for_status();path.write_bytes(response.content)
        assert len(path.read_bytes())==item['size']
        assert item['checksum']=='md5:'+hashlib.md5(path.read_bytes()).hexdigest()
        records.append(dict(file=item['key'],url=item['links']['self'],sha256=sha(path),
            publisher_checksum=item['checksum'],bytes=item['size'],local_cache=path.relative_to(ROOT).as_posix()))
    record={'doi':DOI,'record_id':17200473,'files':records,'executed_external_code':False,
        'license':meta.get('metadata',{}).get('license'),
        'note':'Original downloads stay in ignored .tmp; rerun this script to fetch and verify them.'}
    target=HERE/'reference-source.json'
    if target.exists(): assert data(target)==record, 'Reference bytes changed'
    else:target.write_bytes(encode(record))
    return record

def wl(text):
    return parse_mathematica(text.replace('\\[Mu]','mu')).replace(s.Function('Zeta'),s.zeta).replace(s.Function('PolyLog'),s.polylog)

def raw_W(): return wl((CACHE/'WRen.m').read_text(encoding='utf8'))

def main():
    source=fetch();rows=[]
    def check(name,ok):
        rows.append({'name':name,'passed':bool(ok)});assert ok,name
    def zero(name,x):check(name,s.expand(x)==0)
    result=data(DEST/'result-r001.json')
    r,L,B,CA,CF,nf,TR,a=s.symbols('rho L B CA CF nf TR as4pi')
    ss,tt,mu=s.symbols('s t mu');z,U=s.symbols('z U',positive=True)
    Li2,Li3,Q=s.symbols('Li2 Li3 Q')
    raw=raw_W();x=raw.subs({ss:s.sqrt(z),tt:1,mu:s.sqrt(U),TR:s.Rational(1,2)})
    logs={z:r,1/z:-r,1+z:B,1+1/z:B-r,-z:r+s.I*s.pi,-1/z:-r+s.I*s.pi,
        z*U:L+r/2,U:L-r/2,(z+1)/(z**2*U):B-L-3*r/2,(z+1)/U:B-L+r/2,
        z*U/(z+1):L+r/2-B,U/(z+1):L-r/2-B,(z+1)*U:B+L-r/2}
    nlog=0
    for f in sorted(x.atoms(s.log),key=str):
        matches=[v for k,v in logs.items() if s.cancel(f.args[0]-k)==0]
        check('unique_log_mapping:'+str(f),len(matches)==1)
        x=x.xreplace({f:matches[0]});nlog+=1
    sum3=-Li3-r**2*B/2+r**3/6-s.pi**2*r/3+s.pi**2*B/2+s.I*s.pi*(r*B-B**2-r**2/2)+s.zeta(3)
    polys={(2,-z):Li2,(2,-1/z):-Li2-s.pi**2/6-r**2/2,
        (2,z/(1+z)):-Li2-B**2/2,(2,1/(1+z)):s.pi**2/6+Li2+r*B-B**2/2,
        (3,-z):Li3,(3,-1/z):Li3+(r**3+s.pi**2*r)/6,(3,1+z):Q,(3,1+1/z):sum3-Q}
    npoly=0
    for f in sorted(x.atoms(s.polylog),key=str):
        n,arg=f.args;matches=[v for (nn,k),v in polys.items() if n==nn and s.cancel(arg-k)==0]
        check('unique_polylog_mapping:'+str(f),len(matches)==1)
        x=x.xreplace({f:matches[0]});npoly+=1
    # Verify the differential identities in rho=log(z); constants fixed as in the derivation.
    def dr(f):return s.diff(f,r)+z*s.diff(f,z)+z/(1+z)*s.diff(f,B)-B*s.diff(f,Li2)+Li2*s.diff(f,Li3)
    def ratzero(name,f):check(name,s.cancel(f)==0)
    ratzero('dilog_inversion_derivative',dr(polys[(2,-1/z)])-(B-r))
    ratzero('dilog_fraction_derivative',dr(polys[(2,z/(1+z))])-B/(1+z))
    ratzero('dilog_complement_derivative',dr(polys[(2,1/(1+z))])-z/(1+z)*(r-B))
    ratzero('trilog_inversion_derivative',dr(polys[(3,-1/z)])+polys[(2,-1/z)])
    dsum=(z*(s.pi**2/6-Li2-B*(r+s.I*s.pi))-
          (s.pi**2/6-polys[(2,-1/z)]-(B-r)*(-r+s.I*s.pi)))/(1+z)
    ratzero('trilog_above_cut_sum_derivative',dr(sum3)-dsum)
    x=s.expand(x);W1=s.expand(x.coeff(a,1)/4);W2=s.expand(x.coeff(a,2)/16)
    F2=s.expand(W2-W1**2/2)
    zero('paper_CF_squared_exponentiates',F2.coeff(CF,2))
    check('paper_cut_imaginaries_cancel',not F2.has(s.I))
    check('paper_unresolved_Li3_above_cut_cancels',not F2.has(Q))
    M11=2*s.zeta(3)+(r*r+s.pi**2)*B-s.pi**2*r/2+2*r*Li2-2*Li3-r**3/6
    ours1=CF*s.sympify(result['F1_over_CF']['sympy'])
    ours2=CF*s.sympify(result['F2_over_CF']['sympy']).subs({s.Symbol('T'):nf/2,s.Symbol('M11'):M11})
    zero('paper_one_loop_calibrates_log_dictionary',W1-ours1)
    zero('paper_two_loop_full_finite_function',F2-ours2)
    zero('paper_CA_finite_function',s.expand(F2-ours2).coeff(CA))
    zero('paper_nf_finite_function',s.expand(F2-ours2).coeff(nf))
    # Extract exactly the lambda^2 coefficient of the independent threshold AD file.
    ad=(CACHE/'ADs.m').read_text(encoding='utf8')
    thr=ad.split('\\[Gamma]Thr=',1)[1].split('/.as->4*as4pi;',1)[0]
    thr2=s.expand(wl(re.sub(r'\bas\b','coupling',thr))).coeff(s.Symbol('coupling'),2)
    oursH=s.sympify(result['h2']['sympy']).subs(s.Symbol('T'),nf/2)
    zero('h2_matches_published_threshold_AD',oursH-thr2)
    check('negative_control_compare_W2_to_F2',s.expand(W2-ours2)!=0)
    check('negative_control_wrong_coupling_normalization',s.expand(x.coeff(a,1)-ours1)!=0)
    target=HERE/'reference-checks.json'
    target.write_bytes(encode({'source':source,'checks':rows,'passed':len(rows),
        'dictionary':{'as4pi':'lambda/4','TR':'1/2','nf':'2*T','log(s_paper^2*mu_paper^2)':'L+rho/2',
            'log(t_paper^2*mu_paper^2)':'L-rho/2','s_paper^2/t_paper^2':'z',
            'object':'WRen.m is W_R; subtract W1^2/2 to compare F2=log W_R coefficient'},
        'mapped_logs':nlog,'mapped_polylogs':npoly,'trilog_cut_identity':obj(sum3),
        'analytic_residual':obj(s.expand(F2-ours2)),
        'h2_reference':obj(thr2),'result_sha256':sha(DEST/'result-r001.json'),'code_sha256':sha(__file__),
        'limits':'Positive z, principal polylog boundary values; no general physical-region continuation.'}))
    print(json.dumps({'reference_checks':len(rows),'finite_analytic_residual':0,'h2_residual':0}))

if __name__=='__main__':main()
