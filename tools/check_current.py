"""Standard-library publication/current-contract audit, without recomputing I."""
from pathlib import Path
import argparse,datetime,hashlib,json,re,sys
from urllib.parse import unquote
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();out=a.output_dir.resolve()
    if out.exists():raise ValueError('Use a new output directory')
    if out.is_relative_to(ROOT) and not out.is_relative_to(ROOT/'replay-output'):raise ValueError('Inside the checkout use replay-output/<new>')
    manifest=read(ROOT/'CURRENT-PUBLICATION.json')
    for rel,digest in manifest['files'].items():
        path=(ROOT/rel).resolve();assert path.is_relative_to(ROOT) and sha(path)==digest,rel
    for rel,digest in manifest['external_bindings'].items():assert sha(ROOT/rel)==digest,rel
    from workflow.current_calculation import check_current
    contract=check_current();assert contract['orbits']==163
    selector=read(ROOT/'benchmarks/renormalization/three-loop/current.json');assert selector['full_w3_complete'] and len(selector['channels'])==4
    for row in selector['channels'].values():
        p=(ROOT/'benchmarks/renormalization/three-loop'/row['result']).resolve();assert sha(p)==row['sha256']
        d=read(p)
        for token in row['json_pointer'].strip('/').split('/'):d=d[token]
        assert set(d['U'])==set(map(str,range(-6,1)))
    links=0
    for rel in manifest['active_markdown']+manifest['mutable_navigation']:
        p=ROOT/rel
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8-sig')):
            if '://' in target or target.startswith(('#','mailto:')):continue
            target=unquote(target.strip('<>').split('#',1)[0])
            if target:
                path=(p.parent/target).resolve();assert path.is_relative_to(ROOT) and path.exists(),(rel,target);links+=1
    secret=re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----')
    for rel in manifest['files']:
        p=ROOT/rel;assert p.suffix.lower() not in ['.exe','.dll','.pyc','.pem','.key','.p12','.pfx'],rel
        assert secret.search(p.read_text(encoding='utf-8-sig')) is None,('Credential-like content',rel)
    report={'status':'passed','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'publication_files':len(manifest['files']),'active_local_links':links,'current_contract':contract,'four_color_selector_complete':True,'manifest_sha256':sha(ROOT/'CURRENT-PUBLICATION.json'),'scope':'Publication byte bindings, current runtime contracts, source/result selector and active links. Actual mathematical copied-source replay recorded separately.','formal':'not_run'}
    out.mkdir(parents=True);(out/'checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if not isinstance(v,(dict,list))}))
if __name__=='__main__':main()
