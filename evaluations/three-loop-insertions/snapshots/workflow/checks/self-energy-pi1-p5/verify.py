from pathlib import Path
import json,hashlib,re
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=load(BASE/'prepared/manifest.json')
for d in manifest['inputs']+manifest['outputs']+[manifest['normalization'],manifest['preparer'],manifest['color_backend']]:assert sha(ROOT/d['path'])==d['sha256']
cat=load(BASE/'families/catalog-r007.json');prev=cat['previous_catalog'];assert sha(ROOT/prev['path'])==prev['sha256']
assert sum(o['status']=='computed' for o in cat['orbits'])==13
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
assert (positions,nonzero)==(91,50)
a,b=[load(BASE/f'results/4e-Pi1-{i}/r001.json') for i in ['07','08']]
assert a['I']==b['I'] and a['laurent']==b['laurent']
for file,count in [('exact-checks.json',24),('numeric-checks.json',9)]:
    r=load(HERE/file);assert r['count']==count and all(c['passed'] for c in r['checks'])
rendered=0
for name,doc in [('solutions',BASE/'families/p5-opposite-same/solutions.md'),('07',BASE/'workspaces/4e-Pi1-07/derivation.md'),('08',BASE/'workspaces/4e-Pi1-08/derivation.md')]:
    r=load(HERE/f'render-{name}.json');assert not r['errors'] and r['errorElementsInHtml']==0 and r['sourceSha256']==sha(doc)
    rendered+=r['mathExpressions']
links=0
for p in [BASE/'README.md',BASE/'solve-plan.md',HERE/'README.md',*sorted((BASE/'families').rglob('*.md')),*sorted((BASE/'workspaces').rglob('*.md'))]:
    for m in re.finditer(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
        link=m.group(1)
        if link.startswith(('https://','http://','#')):continue
        assert (p.parent/link.split('#')[0]).resolve().exists(),(p,link)
        links+=1
paths=[HERE/n for n in ['compute.py','numeric.py','publish.py','verify.py','exact-checks.json','numeric-checks.json','render-solutions.json','render-07.json','render-08.json']]
paths+=[BASE/'families/catalog-r007.json',BASE/'families/p5-opposite-same/solutions.md',BASE/'workspaces/4e-Pi1-07/derivation.md',BASE/'workspaces/4e-Pi1-08/derivation.md']
report={'solved':13,'total':13,'positions':positions,'nonzero':nonzero,'exact_checks':24,'numeric_controls':9,'rendered_formulas':rendered,'local_links':links,'results':results,'evidence':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in paths],'formal':'not_run'}
(HERE/'delivery.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(f'13/13 completed; 91 slots / 50 nonzero; {links} links; {rendered} formulas; all hashes verified')
