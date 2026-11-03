"""Check exact bindings, links and preservation, separately from prose review."""
from pathlib import Path
from urllib.parse import unquote
import hashlib,json,re
from prepare import HERE,ROOT,SOURCE,RUNTIME,FAMILY,PUBLIC,read,put,sha,git

def pointer(value,path):
    for part in path.lstrip('/').split('/') if path else []:
        key=part.replace('~1','/').replace('~0','~');value=value[int(key)] if isinstance(value,list) else value[key]
    return value
def walk(value,base):
    total=0
    if isinstance(value,dict):
        if {'path','sha256'}<=value.keys():
            p=(base/value['path']).resolve();assert p.is_file(),p;assert sha(p)==value['sha256'],p
            if 'json_pointer' in value:
                obj=pointer(read(p),value['json_pointer'])
                if 'field_sha256' in value:assert hashlib.sha256(json.dumps(obj,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()==value['field_sha256'],p
            total+=1
        for v in value.values():total+=walk(v,base)
    elif isinstance(value,list):
        for v in value:total+=walk(v,base)
    return total
def main():
    references=0;entries=read(HERE/'entries.json');generated=read(HERE/'evidence/primary/generation/generated.json');cache_entries=0
    for row in entries:
        own=ROOT/row['path'];references+=walk(read(own/'input.json'),own)
        data=read(own/'data/index.json');references+=walk(data,(own/'data'/data['runtime_base']).resolve())
        archived=read(RUNTIME/f"integrals/three-loop/{row['orbit_id']}/subintegrals.json")
        assert read(own/'input.json')['original_ordered_J']==archived['root_J']
        assert read(own/'input.json')['normalization']==archived['normalization']
        assert (own/'data'/data['original_ordered_root']).resolve()==own/'input.json'
        if row['family']=='single-opposite':
            needed={(parameters,str(q)) for parameters,orders in generated[row['orbit_id']]['provenance']['moments'].items() for q in orders}
            indexed=set()
            for unit in data['primitive_units']:
                request=read(RUNTIME/unit['request']['path']);parameters=','.join(map(str,request['parameters']))
                indexed.update((parameters,str(q)) for q in unit['orders'])
            assert indexed==needed,('shared cache dependency mismatch',row['orbit_id'])
            cache_entries+=any(unit['from_shared_batch_cache'] for unit in data['primitive_units'])
        references+=walk(read(own/'checks/index.json'),own/'checks')
        assert sha(own/'result.json')==row['result_sha256']
        assert read(own/'input.json')['normalization'] is not None
        for name in ('README.md','derivation.md','reproduce.md','validation.md','code/reproduce.py','data/index.json','checks/index.json'):assert (own/name).is_file(),(row['orbit_id'],name)
    source=read(SOURCE/'manifest.json')
    for rel,row in source['files'].items():assert sha(SOURCE/rel)==sha(RUNTIME/rel)==row['sha256'],rel
    assert sha(SOURCE/'manifest.json')==sha(RUNTIME/'manifest.json')
    baseline=read(HERE/'preflight.json')
    for rel,h in baseline['protected'].items():assert sha(ROOT/rel)==h,('original drift',rel)
    for rel,h in baseline['public_files'].items():assert sha(PUBLIC/rel)==h,('public drift',rel)
    assert git('rev-parse','HEAD')==baseline['HEAD'] and git('diff','--cached','--name-only').splitlines()==baseline['staged']
    assert git('rev-parse','HEAD',cwd=PUBLIC)==baseline['public_HEAD']
    review=read(HERE/'review/audit.json');assert len(review['entries'])==39
    reviewed=0
    for entry in review['entries']:
        for e in entry['evidence']:
            p=ROOT/e['path'];assert sha(p)==e['file_sha256'],('review drift',p)
            assert e['passage'] in p.read_text(encoding='utf-8-sig').splitlines()[e['line']-1],('passage drift',p)
            reviewed+=1
    docs=[*FAMILY.glob('*.md'),*FAMILY.glob('appendices/**/*.md'),*HERE.rglob('*.md'),ROOT/'integrals/README.md',ROOT/'integrals/three-loop/README.md']
    docs=[p for p in docs if 'runs' not in p.parts and not any(part.startswith('isolation-') for part in p.parts)]
    docs.extend(p for entry in entries for p in (ROOT/entry['path']).rglob('*.md'))
    links=0;broken=[]
    for p in sorted(set(docs)):
        source=p.read_text(encoding='utf-8-sig');source=re.sub(r'\$\$[\s\S]*?\$\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)|(?<!\$)\$(?!\$)[^$\n]*\$', '',source)
        for m in re.finditer(r'(?<!!)\[[^\]\n]+\]\(([^)\n]+)\)',source):
            url=m.group(1).strip().strip('<>')
            if '://' in url or url.startswith(('#','mailto:')):continue
            target=(p.parent/unquote(url.split('#')[0])).resolve();links+=1
            if not target.exists() and target!=HERE/'artifact-check.json':broken.append({'file':p.relative_to(ROOT).as_posix(),'link':url})
    assert not broken,broken
    result={'status':'passed','entries':39,'accepted_results_unchanged':39,'source_and_object_references':references,'single_entries_with_explicit_transitive_cache_dependencies':cache_entries,'runtime_source_files_unchanged':len(read(SOURCE/'manifest.json')['files'])+1,'protected_original_files_unchanged':len(baseline['protected']),'public_files_unchanged':len(baseline['public_files']),'reviewed_passage_bindings_unchanged':reviewed,'current_markdown_files':len(set(docs)),'local_links_checked':links,'broken_links':broken,'HEAD_and_staged_unchanged':True,'scope':'Mechanical source/object bindings and preservation only. Scientific continuity is assessed per entry by the independent reading audit.'}
    put(HERE/'artifact-check.json',result);print(result)
if __name__=='__main__':main()
