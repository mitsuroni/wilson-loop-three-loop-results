"""Verify immutable inputs/results, delivery statuses, links and evidence hashes."""
from pathlib import Path
import json, hashlib, re
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=load(BASE/'prepared/manifest.json')
for item in manifest['inputs']+manifest['outputs']+[manifest['normalization'],manifest['preparer'],manifest['color_backend']]:
    assert sha(ROOT/item['path'])==item['sha256'],item['path']
catalog=load(BASE/'families/catalog-r003.json')
previous=catalog['previous_catalog'];assert sha(ROOT/previous['path'])==previous['sha256']
assert sum(o['status']=='computed' for o in catalog['orbits'])==6
assert catalog['definition_source']['sha256']==sha(BASE/'prepared/normalization.json')
results=[];positions=nonzero=0
for o in catalog['orbits']:
    if o['status']!='computed':continue
    p=BASE/f"results/{o['orbit_id']}/r001.json";r=load(p)
    for key in ['package','manifest']:assert sha(ROOT/r['input'][key])==r['input'][key+'_sha256']
    for dep in r.get('dependencies',[]):assert sha(ROOT/dep['path'])==dep['sha256']
    assert r['formal']=='not_run' and r['normalization']['Euler_count']==1
    c=r['laurent']['coefficients'];assert set(c)==set(map(str,range(-6,1)))
    positions+=len(c);nonzero+=sum(s.sympify(v['expression'])!=0 for v in c.values())
    results.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)})
assert (positions,nonzero)==(42,27)
for left,right in [('01','02'),('05','06')]:
    a=load(BASE/f'results/4e-Pi1-{left}/r001.json');b=load(BASE/f'results/4e-Pi1-{right}/r001.json')
    assert a['I']==b['I'] and a['laurent']==b['laurent']
exact=load(HERE/'exact-checks.json');num=load(HERE/'numeric-checks.json')
assert exact['count']==69 and all(c['passed'] for c in exact['checks'])
assert num['count']==46 and all(c['passed'] for c in num['checks'])
links=0
for p in [BASE/'README.md',BASE/'solve-plan.md',HERE/'README.md',*sorted((BASE/'families').rglob('*.md')),*sorted((BASE/'workspaces').rglob('*.md'))]:
    for match in re.finditer(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
        url=match.group(1)
        if url.startswith(('http://','https://','#')):continue
        assert (p.parent/url.split('#')[0]).resolve().exists(),(p,url)
        links+=1
sources=[HERE/'compute.py',HERE/'numeric.py',HERE/'publish.py',HERE/'verify.py',HERE/'exact-checks.json',HERE/'numeric-checks.json',BASE/'families/catalog-r003.json']
sources+=sorted(HERE.glob('render-*.json'))
for name in ['p1-adjacent-double','p2-adjacent-middle']:sources.append(BASE/f'families/{name}/solutions.md')
for i in ['01','02','05','06']:sources.append(BASE/f'workspaces/4e-Pi1-{i}/derivation.md')
report={'prepared_orbits':13,'solved_orbits':6,'new_results':4,'new_positions':28,'new_nonzero':20,'total_positions':positions,'total_nonzero':nonzero,'exact_checks':69,'numerical_controls':46,'links':links,'results':results,'evidence':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in sources],'formal':'not_run'}
(HERE/'delivery.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(f'6/13 solved; 42 slots / 27 nonzero; {links} links and all source/result hashes verified')
