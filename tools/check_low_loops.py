"""Audit version identity and active navigation of the complete low-loop packet.

Standard library only. This deliberately does not equate hashes with mathematics;
run the separate calculation entrypoints for the exact and numerical checks.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,re
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def within(base,rel):
    p=(base/rel).resolve()
    assert p.is_relative_to(base.resolve()),('Path escapes the declared root',rel)
    return p

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    out=a.output_dir.resolve()
    if out.is_relative_to(ROOT) and not out.is_relative_to(ROOT/'replay-output'):
        raise SystemExit('Inside the checkout, use replay-output/<new directory>.')
    if out.exists():raise SystemExit('Use a new output directory.')
    manifest=read(ROOT/'LOW-LOOP-PUBLICATION.json');count=0
    for rel,digest in manifest['files'].items():
        p=within(ROOT,rel);assert sha(p)==digest,('Changed publication file',rel);count+=1
    for rel,digest in manifest['external_bindings'].items():
        assert sha(within(ROOT,rel))==digest,('Changed upstream binding',rel)
    dm=read(ROOT/'data/manifest.json')
    for rel,digest in dm['files'].items():assert sha(within(ROOT/'data',rel))==digest,('Changed raw/prepared data',rel)
    packet=ROOT/'results/2026-11-02-low-loops';flow=read(packet/'data-flow.json')['integrals']
    assert len(flow)==14 and len({(x['loop_order'],x['orbit_id']) for x in flow})==14
    assert [sum(x['loop_order']==n for x in flow) for n in [1,2]]==[2,12]
    positions=0
    for row in flow:
        assert sha(within(ROOT,row['input']['path']))==row['input']['sha256']
        assert sha(within(ROOT,row['result']))==row['result_sha256']
        assert within(ROOT,row['derivation']).is_file()
        value=read(within(ROOT,row['result']))
        cs=value['laurent']['coefficients'] if row['loop_order']==1 else value['coefficients']
        assert set(map(int,cs))==set(range(row['coefficient_range'][0],row['coefficient_range'][1]+1))
        positions+=len(cs)
    catalogue=read(ROOT/'integrals/two-loop/catalog.json')
    for row in catalogue['integrals']:
        assert row['assembly_sign']==1
        if row['family']=='Pi1':assert row['legacy_to_canonical_multiplier']==-1
        p=within(ROOT/'integrals/two-loop',row['canonical_source']['path'])
        assert sha(p)==row['canonical_source']['sha256']
    assert read(packet/'checks.json')['status']=='passed'
    assert read(packet/'checks.json')['checks']==68
    assert read(packet/'three-loop-budget.json')['status']=='reassembled_and_matched'
    assert read(packet/'numeric-samples.json')['status']=='passed'
    links=0
    active=[ROOT/r for r in manifest['files'] if r.endswith('.md') and '/sources/' not in r]
    active += [ROOT/r for r in manifest['mutable_navigation']]
    for p in active:
        text=p.read_text(encoding='utf-8-sig')
        assert not any(ord(c)<32 and c not in '\t\r\n' for c in text),p
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if '://' in target or target.startswith(('#','mailto:')):continue
            rel=unquote(target.split('#',1)[0].strip('<>'))
            if rel:
                path=(p.parent/rel).resolve()
                assert path.is_relative_to(ROOT) and path.exists(),(p.relative_to(ROOT).as_posix(),rel)
                links+=1
    secret=re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----')
    for rel in manifest['files']:
        p=ROOT/rel
        assert p.suffix.lower() not in ['.pyc','.exe','.dll','.pem','.key','.pfx'],rel
        assert secret.search(p.read_text(encoding='utf-8-sig')) is None,('Credential-like content',rel)
    result={'status':'passed','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'publication_files':count,
     'preserved_upstream_data_files':len(dm['files']),'integrals':len(flow),'individual_Laurent_positions':positions,
     'active_local_links':links,'manifest_sha256':sha(ROOT/'LOW-LOOP-PUBLICATION.json'),
     'scope':'Exact version identities, source/result connections, declared epsilon ranges, active navigation and recorded validation consistency. Not new integral computation or original-physics proof.',
     'formal':'not_run'}
    out.mkdir(parents=True);(out/'checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':main()
