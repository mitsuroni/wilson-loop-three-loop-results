from pathlib import Path
import json,hashlib,re
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1-chain'
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=load(BASE/'prepared/manifest.json')
for d in manifest['inputs']+manifest['outputs']+[manifest['normalization'],manifest['preparer'],manifest['color_backend']]:assert sha(ROOT/d['path'])==d['sha256'],d['path']
cat=load(BASE/'families/catalog-r002.json');assert cat['solved']==2 and cat['prepared']==2
assert sha(ROOT/cat['previous_catalog']['path'])==cat['previous_catalog']['sha256']
positions=nonzero=color_positions=color_nonzero=0;result_refs=[]
for o in cat['orbits']:
    assert o['status']=='computed';ref=o['result'];p=ROOT/ref['path'];assert sha(p)==ref['sha256']
    r=load(p);result_refs.append(ref)
    for name in ['package','manifest']:assert sha(ROOT/r['input'][name])==r['input'][name+'_sha256']
    for d in r['dependencies']:assert sha(ROOT/d['path'])==d['sha256'],d['path']
    assert r['normalization']['Euler_count']==1 and r['normalization']['a1']=='1' and r['formal']=='not_run'
    table=r['laurent']['coefficients'];assert set(table)==set(map(str,range(-6,1)))
    positions+=len(table);nonzero+=sum(s.sympify(v['expression'])!=0 for v in table.values())
    for v in table.values():assert not s.sympify(v['expression']).has(s.EulerGamma)
    assert len(r['color_channels'])==3
    for ch in r['color_channels']:
        assert set(ch['coefficients'])==set(table)
        color_positions+=len(table);color_nonzero+=sum(s.sympify(v['expression'])!=0 for v in ch['coefficients'].values())
    assert r['laurent']['leading_power']==-o['actual_I_pole_order']
assert (positions,nonzero,color_positions,color_nonzero)==(14,9,42,27)
exact=load(HERE/'exact-checks.json');num=load(HERE/'numeric-checks.json')
assert exact['count']==79 and all(c['passed'] for c in exact['checks'])
assert num['count']==48 and all(c['passed'] for c in num['checks'])
# Previous preparation receipts retain their exact artifacts via the archive mapping.
history=load(HERE/'history-manifest.json')['files'];mapping={r['path']:ROOT/r['archived_path'] for r in history}
for r in history:assert sha(ROOT/r['archived_path'])==r['sha256']
old=load(ROOT/'workflow/checks/self-energy-pi1-chain-preparation/delivery.json')
for artifact in old['artifacts']:
    p=mapping.get(artifact['path'],ROOT/artifact['path'])
    assert sha(p)==artifact['sha256'],artifact['path']
rendered=0
for name,p in [('adjacent',BASE/'families/adjacent/solutions.md'),('opposite',BASE/'families/opposite/solutions.md'),
               ('01',BASE/'workspaces/2e-Pi1+Pi1-01/derivation.md'),('02',BASE/'workspaces/2e-Pi1+Pi1-02/derivation.md'),('normalization',BASE/'normalization-and-colors.md')]:
    r=load(HERE/f'render-{name}.json');assert not r['errors'] and r['errorElementsInHtml']==0 and r['sourceSha256']==sha(p)
    rendered+=r['mathExpressions']
links=0
paths=[BASE/'README.md',BASE/'solve-plan.md',BASE/'normalization-and-colors.md',BASE/'prepared/README.md',HERE/'README.md',
       *sorted((BASE/'families').rglob('*.md')),*sorted((BASE/'workspaces').rglob('*.md'))]
for p in paths:
    for m in re.finditer(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
        url=m.group(1)
        if url.startswith(('http:','https:','#')):continue
        assert (p.parent/url.split('#')[0]).resolve().exists(),(p,url);links+=1
files=[HERE/n for n in ['compute.py','numeric.py','publish.py','verify.py','exact-checks.json','numeric-checks.json','first-remainder-range.json','history-manifest.json']]
files+=list(HERE.glob('render-*.json'))+[BASE/'families/catalog-r002.json']+paths
report={'status':'computationally_complete','solved':2,'I_positions':14,'I_nonzero':9,'color_positions':42,'color_nonzero':27,
        'exact_checks':79,'numeric_controls_or_moments':48,'local_links':links,'rendered_expressions':rendered,
        'source_and_prepared_unchanged':True,'previous_preparation_artifacts_preserved':True,'formal':'not_run','results':result_refs,
        'artifacts':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in files]}
(HERE/'delivery.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(f'2/2 I completed; 14 slots/9 nonzero, color 42/27; 79 exact, 48 numeric; {links} links, {rendered} formulas; hashes/history verified')
