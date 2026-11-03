"""Bind initial graph/C/J fields and original ordered domains in each reader input."""
from pathlib import Path
from prepare import ROOT,PUBLIC,RUNTIME,READING,HERE,read,put,sha,exact
def ptr(obj,p):
    for token in p.strip('/').split('/') if p else []:
        token=token.replace('~1','/').replace('~0','~');obj=obj[int(token)] if isinstance(obj,list) else obj[token]
    return obj
def main():
    current={x['orbit_id']:x for x in read(PUBLIC/'benchmarks/current/catalog.json')['entries']};updated=[]
    for row in read(HERE/'entries.json'):
        oid=row['orbit_id'];loop=row['loop_directory'];num=1 if loop=='one-loop' else 2
        if loop=='one-loop':continue # its original input already binds graph/native/ordered/prepared.
        up=read(PUBLIC/f'data/2026-11-02/integrals/{num}loop/{oid}/input.json')
        refs={}
        def walk(v):
            if isinstance(v,dict):
                if 'path' in v and 'sha256' in v and isinstance(v['path'],str):refs[v['path']]=v['sha256']
                for x in v.values():walk(x)
            elif isinstance(v,list):
                for x in v:walk(x)
        walk(up)
        for rel,h in refs.items():
            source=PUBLIC/'data/2026-11-02'/rel
            assert sha(source)==h;exact(source,RUNTIME/'data/2026-11-02'/rel)
        base='../../families/low-loop-runtime/'
        def source_ref(ref):
            return {**ref,'path':base+'data/2026-11-02/'+ref['path']}
        own=READING/loop/oid;input=read(own/'input.json');entry=current[oid];native_ref=entry['native_prepared']
        exact(PUBLIC/native_ref['path'],RUNTIME/native_ref['path']);assert sha(RUNTIME/native_ref['path'])==native_ref['sha256']
        native=ptr(read(RUNTIME/native_ref['path']),native_ref['field'])
        ordered=native.get('J_raw_s',native['J']);order=ordered.get('iterated_order_outer_to_inner',ordered.get('integration_order_outer_to_inner',ordered.get('variables')))
        input.update({'graph_definition':source_ref(up['graph_definition']),
            'native_raw_orbit':source_ref(up['native_raw_orbit']),
            'native_prepared':{**native_ref,'path':base+native_ref['path']},
            'original_ordered_J':{'path':base+native_ref['path'],'json_pointer':native_ref['field']+('/J_raw_s' if 'J_raw_s' in native else '/J'),'sha256':native_ref['sha256'],
                'order_outer_to_inner':order,'bounds':ordered.get('bounds'),'alpha_measure':ordered.get('alpha_measure'),'jacobian':ordered.get('jacobian')},
            'normalization':native['normalization'],'operator':native.get('operator'),
            'scope':'Initial graph/raw/native C/J and original ordered domain bound explicitly; accepted normalized I and historical replay comparison remain separate.'})
        put(own/'input.json',input);updated.append(oid)
    files={p.relative_to(RUNTIME).as_posix():sha(p) for p in RUNTIME.rglob('*') if p.is_file() and p!=RUNTIME/'source-manifest.json'}
    m=read(RUNTIME/'source-manifest.json');m['files']=files;m['initial_source_binding_update']={'orbits':updated,'scope':'Additional actual initial-source objects; no original runtime object modified.'};put(RUNTIME/'source-manifest.json',m)
    put(HERE/'source-binding-check.json',{'status':'passed','new_two_loop_initial_bindings':len(updated),'orbits':updated,'runtime_files':len(files)})
    print({'initial_bindings_added':len(updated),'runtime_files':len(files)})
if __name__=='__main__':main()
