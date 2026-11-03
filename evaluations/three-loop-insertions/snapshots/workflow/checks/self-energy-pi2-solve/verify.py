from pathlib import Path
import json,hashlib,re
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def verify_records(records):
    for entry in records:assert sha(ROOT/entry['path'])==entry['sha256'],entry['path']
def main():
    prep=load(HERE/'preparation-checks.json');exact=load(HERE/'exact-checks.json');numeric=load(HERE/'numeric-checks.json')
    assert (prep['count'],exact['count'],numeric['count'])==(50,166,33)
    for report in [prep,exact,numeric]:assert all(c['passed'] for c in report['checks'])
    assert exact['positions']==42 and exact['nonzero']==25
    tests=load(HERE/'tests-status.json');assert tests['exit_code']==0
    meta=load(BASE/'prepared/manifest.json');verify_records(meta['inputs']+meta['outputs']+meta['code'])
    verify_records(numeric['dependencies'])
    for oid in ['2e-Pi2-01','2e-Pi2-02']:
        result=load(BASE/f'results/{oid}/r001.json');verify_records(result['dependencies'])
        assert sha(ROOT/result['input']['package'])==result['input']['package_sha256']
        assert sha(ROOT/result['input']['manifest'])==result['input']['manifest_sha256']
        assert len(result['color_channels'])==3
        for channel in result['color_channels']:
            rows=channel['laurent']['coefficients']
            assert set(rows)==set(map(str,range(-6,1)))
            assert all(not any(token in value['expression'] for token in ['gamma','EulerGamma','eps','p_AA','p_Af','p_Ff']) for value in rows.values())
        files=list((BASE/f'workspaces/{oid}').iterdir())
        assert len(files)==1 and files[0].name=='derivation.md'
    # Preserve the previous production delivery's immutable navigation artifact.
    history=load(HERE/'history-manifest.json')
    old=load(ROOT/'workflow/checks/self-energy-pi2-montopar/delivery.json')
    for entry in old['artifacts']:
        p=HERE/'history'/entry['path'] if entry['path'] in history else ROOT/entry['path']
        assert sha(p)==entry['sha256'],entry['path']
    doc_reports=[('families/adjacent/solutions.md','render-adjacent.json'),
                 ('families/opposite/solutions.md','render-opposite.json'),('operator-expansion.md','render-operator.json')]
    expressions=0
    for doc,report in doc_reports:
        r=load(HERE/report);assert r['sourceSha256']==sha(BASE/doc) and not r['errors'] and r['errorElementsInHtml']==0
        expressions+=r['mathExpressions']
    docs=[BASE/doc for doc,_ in doc_reports]+[BASE/'README.md',BASE/'families/current-results.md',BASE/'prepared/README.md',HERE/'README.md']
    docs+=list((BASE/'workspaces').glob('*/derivation.md'))
    links=0
    for p in docs:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if link.startswith(('https:','http:','#')):continue
            target=(p.parent/link.split('#')[0]).resolve()
            if target==HERE/'delivery.json':continue
            assert target.exists(),(p,link)
            links+=1
    paths=[*docs,HERE/'prepare.py',HERE/'compute.py',HERE/'numeric.py',HERE/'publish.py',HERE/'test_prepare.py',HERE/'verify.py',HERE/'preparation-checks.json',HERE/'exact-checks.json',HERE/'numeric-checks.json',HERE/'tests-status.json',BASE/'operator-expansion-r001.json',BASE/'families/catalog-r001.json']
    report={'status':'computed_through_finite','orbits':2,'color_channels':6,'positions':42,'nonzero':25,
        'preparation_checks':50,'exact_checks':166,'numeric_controls':33,'unittest_methods':3,'invalid_input_subcases':7,
        'rendered_expressions':expressions,'local_links':links,'history_verified':True,'formal':'not_run',
        'artifacts':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in paths]}
    (HERE/'delivery.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='artifacts'},ensure_ascii=False))
if __name__=='__main__':main()
