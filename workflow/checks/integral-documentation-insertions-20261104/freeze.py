"""Inventory the final new reader, exact runtime and executed evidence."""
from prepare import HERE,ROOT,FAMILY,read,put,sha
def main():
    files={p for p in FAMILY.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    for e in read(HERE/'entries.json'):files.update(p for p in (ROOT/e['path']).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in HERE.rglob('*') if p.is_file() and p!=HERE/'manifest.json' and 'runs' not in p.parts and not any(v.startswith('isolation-r') for v in p.relative_to(HERE).parts[:-1]) and '__pycache__' not in p.parts)
    # Other active research windows also update the shared roadmap. Bind our
    # paragraph separately, instead of claiming its complete mutable file frozen.
    navigation=[]
    for rel in ('workflow/master-plan.md','workflow/project-map.md','workflow/data-flow.md'):
        paragraphs=(ROOT/rel).read_text(encoding='utf-8-sig').split('\n\n')
        matching=[p for p in paragraphs if ('WILSON-5.27' in p or '15个二点／交换插入当前正文' in p)]
        assert len(matching)==1,(rel,len(matching));navigation.append({'path':rel,'task_paragraph':matching[0]})
    put(HERE/'navigation-update.json',{'scope':'Exact task-owned paragraphs in shared mutable navigation; unrelated concurrent four-leg progress is retained. Complete navigation files are not frozen by this manifest.','paragraphs':navigation});files.add(HERE/'navigation-update.json')
    for rel in ('integrals/README.md','integrals/three-loop/README.md','workflow/backlog/tasks/wilson-5.27 - Complete-continuous-derivations-and-audit-fifteen-three-loop-insertion-integrals.md'):files.add(ROOT/rel)
    put(HERE/'manifest.json',{'format':'wilson-insertion-continuous-reader/1','entries':15,'scope':'New reader, exact old-runtime copy, source/object indices and actual evidence; temporary copies/caches excluded. This inventory grants no mathematical proof.','files':{p.relative_to(ROOT).as_posix():{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(files)}})
    print({'manifest_files':len(files)})
if __name__=='__main__':main()
