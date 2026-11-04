"""Mechanical admission of exact bindings and links; prose review is separate."""
from urllib.parse import unquote
import json,hashlib,re
from prepare import HERE,ROOT,SOURCE,FAMILY,RUNTIME,PUBLIC,read,put,sha,ptr,git
def walk(value,base):
    count=0
    if isinstance(value,dict):
        if {'path','sha256'}<=value.keys():
            p=(base/value['path']).resolve();assert p.is_file(),p;assert sha(p)==value['sha256'],p
            pointer=value.get('json_pointer',value.get('field'))
            if pointer is not None:
                field=ptr(read(p),pointer)
                if 'field_sha256' in value:assert hashlib.sha256(json.dumps(field,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()==value['field_sha256'],p
            count+=1
        for v in value.values():count+=walk(v,base)
    elif isinstance(value,list):
        for v in value:count+=walk(v,base)
    return count
def main():
    entries=read(HERE/'entries.json');refs=0
    for e in entries:
        own=ROOT/e['path'];i=read(own/'input.json');d=read(own/'data/index.json')
        refs+=walk(i,own)+walk(d,RUNTIME)+walk(read(own/'checks/index.json'),own/'checks')
        assert sha(own/'result.json')==e['result_sha256'];assert i['original_ordered_J'] is not None
        native=ptr(read(own/i['source_bindings']['native_prepared']['path']),i['source_bindings']['native_prepared']['field'])
        assert i['native_normalization']==native['normalization']
        old=read(RUNTIME/f"integrals/three-loop/{e['orbit_id']}/subintegrals.json")
        expected=old.get('raw_ordered_J',native.get('J_raw_s',native.get('J')))
        assert i['original_ordered_J']==expected and i['evaluated_quantity_normalization']==read(own/'result.json')['normalization']
        assert (own/'data'/d['original_ordered_root']).resolve()==own/'input.json'
        for name in ('README.md','derivation.md','reproduce.md','validation.md','code/reproduce.py','data/index.json','checks/index.json'):assert (own/name).is_file(),(e['orbit_id'],name)
    manifest=read(SOURCE/'manifest.json')
    for rel,row in manifest['files'].items():assert sha(SOURCE/rel)==sha(RUNTIME/rel)==row['sha256'],rel
    assert sha(SOURCE/'manifest.json')==sha(RUNTIME/'manifest.json')
    base=read(HERE/'preflight.json')
    for rel,h in base['protected'].items():assert sha(ROOT/rel)==h,('original drift',rel)
    for rel,h in base['public_files'].items():assert sha(PUBLIC/rel)==h,('public drift',rel)
    assert git('rev-parse','HEAD')==base['HEAD'] and git('diff','--cached','--name-only').splitlines()==base['staged']
    assert git('rev-parse','HEAD',cwd=PUBLIC)==base['public_HEAD']
    review=read(HERE/'review/audit.json');assert len(review['entries'])==15;reviewed=0
    for e in review['entries']:
        for ev in e['evidence']:
            p=ROOT/ev['path'];assert sha(p)==ev['file_sha256'],('review drift',p)
            assert ev['passage'] in p.read_text(encoding='utf-8-sig').splitlines()[ev['line']-1],('passage drift',p);reviewed+=1
    docs=[*FAMILY.glob('*.md'),*FAMILY.glob('appendices/**/*.md'),*HERE.rglob('*.md'),ROOT/'integrals/README.md',ROOT/'integrals/three-loop/README.md']
    docs=[p for p in docs if 'runs' not in p.parts and not any(part.startswith('isolation-r') for part in p.parts)]
    docs.extend(p for e in entries for p in (ROOT/e['path']).rglob('*.md'))
    broken=[];links=0
    for p in sorted(set(docs)):
        text=p.read_text(encoding='utf-8-sig');text=re.sub(r'\$\$[\s\S]*?\$\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)|(?<!\$)\$(?!\$)[^$\n]*\$', '',text)
        for m in re.finditer(r'(?<!!)\[[^\]\n]+\]\(([^)\n]+)\)',text):
            url=m.group(1).strip().strip('<>')
            if '://' in url or url.startswith(('#','mailto:')):continue
            target=(p.parent/unquote(url.split('#')[0])).resolve();links+=1
            if not target.exists() and target!=HERE/'artifact-check.json':broken.append({'file':p.relative_to(ROOT).as_posix(),'link':url})
    assert not broken,broken
    result={'status':'passed','entries':15,'adopted_results_unchanged':15,'exact_source_and_object_references':refs,'runtime_files_unchanged':len(manifest['files'])+1,'protected_original_files_unchanged':len(base['protected']),'public_files_unchanged':len(base['public_files']),'reviewed_passage_bindings_unchanged':reviewed,'current_markdown_files':len(set(docs)),'local_links_checked':links,'broken_links':broken,'HEAD_and_staged_unchanged':True,'scope':'Mechanical identities and preservation, not a proof or an automatic scientific continuity grade.'}
    put(HERE/'artifact-check.json',result);print(result)
if __name__=='__main__':main()
