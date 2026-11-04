from pathlib import Path
import json,hashlib,re
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=load(BASE/'prepared/manifest.json')
for d in manifest['inputs']+manifest['outputs']+[manifest['normalization'],manifest['preparer'],manifest['color_backend']]:assert sha(ROOT/d['path'])==d['sha256']
cat=load(BASE/'families/catalog-r004.json');prev=cat['previous_catalog'];assert sha(ROOT/prev['path'])==prev['sha256']
assert sum(o['status']=='computed' for o in cat['orbits'])==7
positions=nonzero=0;results=[]
for o in cat['orbits']:
    if o['status']!='computed':continue
    p=BASE/f"results/{o['orbit_id']}/r001.json";r=load(p)
    for k in ['package','manifest']:assert sha(ROOT/r['input'][k])==r['input'][k+'_sha256']
    for d in r.get('dependencies',[]):assert sha(ROOT/d['path'])==d['sha256']
    assert r['formal']=='not_run' and r['normalization']['Euler_count']==1
    c=r['laurent']['coefficients'];assert set(c)==set(map(str,range(-6,1)))
    positions+=len(c);nonzero+=sum(s.sympify(v['expression'])!=0 for v in c.values())
    results.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)})
assert (positions,nonzero)==(49,30)
for file,count in [('exact-checks.json',15),('numeric-checks.json',21)]:
    r=load(HERE/file);assert r['count']==count and all(c['passed'] for c in r['checks'])
rendered=0
for name,doc in [('solutions',BASE/'families/p6-opposite-distinct/solutions.md'),('derivation',BASE/'workspaces/4e-Pi1-11/derivation.md')]:
    r=load(HERE/f'render-{name}.json');assert not r['errors'] and r['errorElementsInHtml']==0 and r['sourceSha256']==sha(doc)
    rendered+=r['mathExpressions']
links=0
for p in [BASE/'README.md',BASE/'solve-plan.md',HERE/'README.md',*sorted((BASE/'families').rglob('*.md')),*sorted((BASE/'workspaces').rglob('*.md'))]:
    for m in re.finditer(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
        link=m.group(1)
        if link.startswith(('https://','http://','#')):continue
        assert (p.parent/link.split('#')[0]).resolve().exists(),(p,link)
        links+=1
paths=[HERE/n for n in ['compute.py','numeric.py','publish.py','verify.py','exact-checks.json','numeric-checks.json','render-solutions.json','render-derivation.json']]
paths+=[BASE/'families/catalog-r004.json',BASE/'families/p6-opposite-distinct/solutions.md',BASE/'workspaces/4e-Pi1-11/derivation.md']
report={'solved':7,'total':13,'positions':positions,'nonzero':nonzero,'exact_checks':15,'numeric_controls':21,'rendered_formulas':rendered,'local_links':links,'results':results,'evidence':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in paths],'formal':'not_run'}
(HERE/'delivery.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(f'7/13 completed; 49 slots / 30 nonzero; {links} links; {rendered} formulas; all input/result/dependency hashes verified')
