"""Audit the published source/data snapshot using only the Python standard library.

This verifies exact bytes, catalogue pointers and active navigation. It does not
evaluate the integrals or establish their original mathematical correctness.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, re, sys
from functools import lru_cache
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

@lru_cache(maxsize=32)
def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def within(base,relative):
    p=(base/relative).resolve()
    if not p.is_relative_to(base.resolve()):
        raise ValueError('Published path escapes its root: '+relative)
    return p

def pointer(value,path):
    if path=='':return value
    if not path.startswith('/'):
        raise ValueError('Invalid JSON pointer: '+path)
    for token in path[1:].split('/'):
        token=token.replace('~1','/').replace('~0','~')
        value=value[int(token)] if isinstance(value,list) else value[token]
    return value

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,required=True,help='New directory for audit output; input snapshot is never changed')
    args=ap.parse_args()
    dest=args.output_dir.resolve()
    if dest.exists():
        raise FileExistsError('Use a new audit output directory: '+str(dest))
    counts={'package_source_files':0,'explicit_source_adaptations':0,'recorded_package_tests':0,'data_manifest_files':0,'data_source_assets':0,'integral_records':0,'file_references':0,'field_references':0,'active_local_links':0}
    checked={}
    def bound(p,expected):
        key=p.resolve()
        actual=checked.setdefault(key,sha(p)) if key not in checked else checked[key]
        assert actual==expected,(p.relative_to(ROOT).as_posix(),expected,actual)
    for name in ['qgraf-pipeline','montopar']:
        base=ROOT/'packages'/name
        d=load_json(str(base/'PUBLICATION-SOURCE.json'))
        for original,row in d['source_files'].items():
            p=within(base,row['published_path'])
            bound(p,row['published_sha256'])
            if row['byte_identical']:
                assert row['published_sha256']==row['original_sha256'],original
            else:
                counts['explicit_source_adaptations']+=1
                assert d.get('publication_adapters'),original
            counts['package_source_files']+=1
        for p in base.rglob('*'):
            if p.is_file():
                assert not any(x in p.parts for x in ['__pycache__','.pytest_cache','.git','.superpowers']),p
                assert p.name not in ['CLAUDE.md','AGENTS.md','qgraf_paths.yaml'],p
                assert p.suffix not in ['.pyc','.pyo','.exe','.dll'],p
        record=load_json(str(base/'VALIDATION.json'))
        assert record['status']=='passed' and record['failures']==record['errors']==0,name
        assert record['source_manifest_sha256']==sha(base/'PUBLICATION-SOURCE.json'),name
        counts['recorded_package_tests']+=record['original_tests_executed']
    version=ROOT/'data/2026-11-02'
    catalog=load_json(str(version/'catalog.json'))
    manifest=load_json(str(ROOT/'data/manifest.json'))
    for relative,digest in manifest['files'].items():
        bound(within(ROOT/'data',relative),digest)
        counts['data_manifest_files']+=1
    for row in catalog['source_assets']:
        bound(within(version,row['path']),row['sha256'])
        counts['data_source_assets']+=1
    def refs(value):
        if isinstance(value,dict):
            if isinstance(value.get('path'),str) and isinstance(value.get('sha256'),str):
                p=within(version,value['path'])
                bound(p,value['sha256'])
                counts['file_references']+=1
                if 'json_pointer' in value:
                    field=pointer(load_json(str(p)),value['json_pointer'])
                    if value.get('field_sha256'):
                        field_sha=hashlib.sha256(json.dumps(field,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode('utf-8')).hexdigest()
                        assert field_sha==value['field_sha256'],(p,value['json_pointer'])
                    counts['field_references']+=1
            for v in value.values():refs(v)
        elif isinstance(value,list):
            for v in value:refs(v)
    refs(catalog)
    ids=[]
    for row in catalog['integrals']:
        ids.append(row.get('orbit_id',row.get('id')))
        inp=within(version,row['input_path'])
        readme=within(version,row['readme_path'])
        assert inp.is_file() and readme.is_file(),row
        refs(load_json(str(inp)))
        counts['integral_records']+=1
    assert len(ids)==len(set(ids))==163
    active=[ROOT/'README.md',ROOT/'packages/README.md',ROOT/'packages/qgraf-pipeline/README.md',ROOT/'packages/montopar/README.md',ROOT/'data/README.md',ROOT/'integrals/README.md',ROOT/'results/README.md',ROOT/'templates/integral-derivation.md',ROOT/'docs/reproduction.md',ROOT/'docs/publication-status.md']
    active+=list((ROOT/'docs/upstream').rglob('*.md'))
    active+=list((ROOT/'validation/upstream-2026-11-02').glob('*.md'))
    if (version/'README.md').exists():active.append(version/'README.md')
    active+=list((version/'integrals').rglob('README.md'))
    missing=[]
    for p in active:
        text=p.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\t\r\n' for c in text),p
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if '://' in target or target.startswith(('#','mailto:')):continue
            target=unquote(target.split('#',1)[0].strip('<>'))
            if target:
                counts['active_local_links']+=1
                if not (p.parent/target).resolve().exists():
                    missing.append([p.relative_to(ROOT).as_posix(),target])
    assert not missing,missing[:30]
    secret_pattern=re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----')
    for folder in ['packages','data']:
        for p in (ROOT/folder).rglob('*'):
            if p.is_file():
                assert p.suffix not in ['.pem','.key','.p12','.pfx'],p
                text=p.read_text(encoding='utf-8-sig')
                assert secret_pattern.search(text) is None,('Credential-like content: file withheld',p.relative_to(ROOT).as_posix())
    result={'status':'passed','run_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'counts':counts,'missing_active_links':missing,'inputs_byte_identical_to_published_manifests':True,'scope':'Source/data integrity, catalogue field bindings and active reader navigation; not a new proof of graph completeness or integral correctness. Tests and representative CLI replays are recorded separately.','formal':'not_run'}
    dest.mkdir(parents=True)
    (dest/'checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':
    main()
