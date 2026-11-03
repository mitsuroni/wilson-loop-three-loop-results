"""Bind current metadata explicitly, without editing archived replay inputs."""
from prepare import HERE,ROOT,READING,read,put

def main():
    for row in read(HERE/'entries.json'):
        own=ROOT/row['path'];obj=read(own/'input.json')
        if row['loop_directory']=='one-loop':
            def bind(value):
                if isinstance(value,dict):
                    if isinstance(value.get('path'),str) and value['path'].startswith(('data/','integrals/')):
                        value['path']='../../families/low-loop-runtime/'+value['path']
                    for v in value.values():bind(v)
                elif isinstance(value,list):
                    for v in value:bind(v)
            bind(obj);obj['references_relative_to']='this input.json directory'
        else:obj['references_relative_to']='this input.json directory; family_catalog_entry is archived catalog content, relative to its original two-loop runtime directory'
        put(own/'input.json',obj)
        checks=read(own/'checks/index.json')
        checks['scientific_document_review']='../../../../workflow/checks/integral-documentation-low-loops-20261103/review/README.md'
        checks['scope']='Actual replay, archived controls and independent scientific document review have separate reports. Replay does not certify prose continuity or rerun all historical CAS integrations.'
        put(own/'checks/index.json',checks)
    print({'current_metadata_bound':14,'archived_inputs_changed':0})
if __name__=='__main__':main()
