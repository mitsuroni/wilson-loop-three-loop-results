"""Mechanical bindings, links and preservation; never a mathematical prose grade."""
from pathlib import Path
import re,json
from urllib.parse import unquote
from prepare import HERE,ROOT,PUBLIC,READING,RUNTIME,read,put,sha,git

def pointer(value,path):
    for part in path.lstrip('/').split('/') if path else []:
        key=part.replace('~1','/').replace('~0','~')
        value=value[int(key)] if isinstance(value,list) else value[key]
    return value

def main():
    refs=0;entry_results=0
    for row in read(HERE/'entries.json'):
        own=ROOT/row['path'];obj=read(own/'input.json')
        for key,ref in obj.items():
            if not isinstance(ref,dict) or not {'path','sha256'}<=ref.keys():continue
            p=(own/ref['path']).resolve();assert p.is_file(),(row['orbit_id'],key,p)
            assert sha(p)==ref['sha256'],(row['orbit_id'],key,'hash')
            field=ref.get('json_pointer',ref.get('field'))
            if field is not None:
                value=pointer(read(p),field)
                if 'field_sha256' in ref:
                    import hashlib
                    h=hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
                    assert h==ref['field_sha256'],(row['orbit_id'],key,'field')
            refs+=1
        assert sha(own/'result.json')==row['result_sha256'];entry_results+=1
        data=read(own/'data/index.json');base=(own/'data'/data['runtime_base']).resolve()
        for category in ('adopted_primitives','assembly_code','comparison_targets'):
            for ref in data[category]:
                p=base/ref['path'];assert p.is_file() and sha(p)==ref['sha256'],(row['orbit_id'],category,ref['path'])
                refs+=1
        for name in ('README.md','derivation.md','reproduce.md','validation.md','code/reproduce.py','data/index.json','checks/index.json'):
            assert (own/name).is_file(),(row['orbit_id'],name)
    runtime=read(RUNTIME/'source-manifest.json')['files']
    for rel,h in runtime.items():assert sha(RUNTIME/rel)==h,rel
    baseline=read(HERE/'preflight.json')
    for rel,h in baseline['protected'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in baseline['public_files'].items():assert sha(PUBLIC/rel)==h,rel
    assert git('rev-parse','HEAD')==baseline['HEAD'];assert git('diff','--cached','--name-only').splitlines()==baseline['staged']
    assert git('rev-parse','HEAD',cwd=PUBLIC)==baseline['public_HEAD']
    skill=ROOT/'workflow/skills/wilson-integral-derivations';installed=Path('C:/Users/Lenovo/.codex/skills/wilson-integral-derivations')
    for rel,h in read(HERE/'skill-install.json')['files'].items():assert sha(skill/rel)==h==sha(installed/rel),rel
    review=read(HERE/'review/audit.json');review_bindings=0
    for entry in review['entries']:
        for evidence in entry['evidence']:
            p=ROOT/evidence['path'];assert sha(p)==evidence['file_sha256'],('reviewed article drift',p)
            assert evidence['passage'] in p.read_text(encoding='utf-8-sig').splitlines()[evidence['line']-1],('reviewed passage drift',p)
            review_bindings+=1
        original=entry['original_source_evidence'];assert sha(ROOT/original['input_path'])==original['input_sha256']
    jet=read(HERE/'review/hlog-transcription-check.json');assert sha(ROOT/jet['markdown'])==jet['markdown_sha256']
    docs=[READING/'README.md',READING/'one-loop/README.md',READING/'two-loop/README.md',RUNTIME/'README.md',ROOT/'workflow/integral-documentation-standard.md',*skill.rglob('*.md'),*HERE.rglob('*.md')]
    docs=[p for p in docs if 'isolation-r001' not in p.parts and 'runs' not in p.parts]
    docs.extend(p for folder in ('one-loop','two-loop') for p in (READING/folder).rglob('*.md'))
    broken=[];links=0
    for p in sorted(set(docs)):
        source=p.read_text(encoding='utf-8-sig')
        # Formula coefficient notation such as [e^j](...) is not a prose link.
        source=re.sub(r'\$\$[\s\S]*?\$\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)|(?<!\$)\$(?!\$)[^$\n]*\$', '',source)
        for match in re.finditer(r'(?<!!)\[[^\]\n]+\]\(([^)\n]+)\)',source):
            url=match.group(1).strip().strip('<>')
            if '://' in url or url.startswith(('#','mailto:')):continue
            url=unquote(url.split('#')[0]);target=(p.parent/url).resolve();links+=1
            if not target.exists() and target!=HERE/'artifact-check.json':broken.append({'file':p.relative_to(ROOT).as_posix(),'link':url})
    assert not broken,broken
    result={'status':'passed','current_entries':14,'accepted_results_unchanged':entry_results,'source_references_checked':refs,'runtime_files_unchanged':len(runtime),
        'protected_original_files_unchanged':len(baseline['protected']),'public_files_unchanged':len(baseline['public_files']),
        'current_markdown_files':len(set(docs)),'local_links_checked':links,'broken_links':broken,'installed_skill_same_bytes':True,'reviewed_passage_bindings_unchanged':review_bindings,'HEAD_and_staged_unchanged':True,
        'scope':'Mechanical existence, exact bindings, links and immutable-source preservation only. Scientific continuity is graded by the independent review, not by this script.'}
    put(HERE/'artifact-check.json',result);print(result)
if __name__=='__main__':main()
