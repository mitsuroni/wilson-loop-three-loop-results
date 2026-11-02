"""Read and replay the single current scientific calculation, never legacy signs.

This is a bound reader/launcher for existing data and algebra engines, not a new
integral IR or proof system. The default47 ordered solver roots retain their role.
"""
from pathlib import Path
import argparse,copy,hashlib,json,re,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
CONVENTION='wilson-qcd-action-path-fourier/1.0'
MANIFEST=ROOT/'workflow/current-calculation.json'
def _json(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def _sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def _path(root,rel):
    relpath=Path(rel)
    if relpath.is_absolute():raise ValueError('Current bindings require repository-relative paths')
    p=(root/relpath).resolve()
    if not p.is_relative_to(root.resolve()):raise ValueError('Current binding escapes repository')
    return p
def _bound(root,ref):
    p=_path(root,ref['path'])
    if _sha(p)!=ref['sha256']:raise ValueError('Changed current input: '+ref['path'])
    return _json(p)
def _header(value):
    if any(k in value for k in ['legacy_sign','family_sign','legacy_to_canonical_multiplier','legacy_mode']):
        raise ValueError('Legacy conversion fields are not part of the current interface')
    if value.get('physical_convention')!=CONVENTION:raise ValueError('Current physical convention required')
    if value.get('family_sign_conversion_allowed') is not False:raise ValueError('Legacy family sign conversion is forbidden')
    if value.get('assembly_sign')!=1:raise ValueError('Extra current assembly sign is forbidden')
def _native(value):
    if value.get('physical_convention')!=CONVENTION:raise ValueError('Native C/J convention required')
    for row in value['orbits']:
        norm=row.get('normalization',{});prep=row.get('preparation',{})
        contract=row.get('physical_convention',{})
        expected={'id':CONVENTION,'wilson_transport':"U'=ig*A*U",'momentum_Fourier':'exp(-ipX)',
                  'momentum_Wick':'I','coordinate_Wick':'-I','Dyson_external_coefficient':1,
                  'raw_Euler_count':0,'vertex_labels_unchanged':True,'one_pi_includes_internal_loop_signs':True}
        if any(contract.get(key)!=v for key,v in expected.items()) or norm.get('physical_convention')!=CONVENTION:
            raise ValueError('Mixed or undeclared orbit-level physical convention')
        chronological=contract.get('chronological_word')
        if not isinstance(chronological,list) or contract.get('wilson_word')!=list(reversed(chronological)):
            raise ValueError('Wilson word does not follow the current transport equation')
        if norm.get('assembly_sign')!=1 or norm.get('euler_compensation_count')!=1:
            raise ValueError('Native C/J sign/Euler normalization mismatch')
        if prep.get('native_sign_conversion_count')!=0 or prep.get('J_changed') is not False:
            raise ValueError('Native preparation must preserve J without sign conversion')
def _pointer(value,field):
    for token in field.strip('/').split('/') if field else []:
        token=token.replace('~1','/').replace('~0','~')
        value=value[int(token)] if isinstance(value,list) else value[token]
    return value

def _entry_definition(row,package):
    if package.get('loop_order')!=row['loop_order'] or package.get('emission_type')!=row['family']:
        raise ValueError('Current orbit loop/family does not match its source package')
    value=_pointer(package,row['native_prepared']['field'])
    if value.get('orbit_id')!=row['orbit_id']:
        raise ValueError('Current catalog field points to another orbit')
    return value

def _whole_i_views(candidates,data):
    byfield={x['field']:x for x in candidates}
    if isinstance(data.get('whole_I_coefficients'),dict):
        powers=[]
        for key in data['whole_I_coefficients']:
            match=re.fullmatch(r'(?:eps\^)?(-?\d+)',key)
            if not match:raise ValueError('Unknown explicit whole-I Laurent key')
            powers.append(int(match[1]))
        byfield['/whole_I_coefficients']={'field':'/whole_I_coefficients','powers':sorted(powers)}
    for field in ['/whole_I_laurent/coefficients','/whole_I_coefficients']:
        if field in byfield:return [byfield[field]|{'quantity_role':'whole normalized I'}]
    channels=[x|{'quantity_role':'whole normalized I','channel':x['field'].split('/')[2]}
              for x in candidates if re.fullmatch(r'/channels/(?:g|q)/I_coefficients',x['field'])]
    if channels:return channels
    if '/I_geom_laurent/coefficients' in byfield:
        return [byfield['/I_geom_laurent/coefficients']|{'quantity_role':'shared scalar I_geom; exact external operator/color still required'}]
    for field in ['/laurent/coefficients','/laurent_coefficients','/laurent_Hlog_coefficients','/laurent/coefficients_input_form']:
        if field in byfield:return [byfield[field]|{'quantity_role':'whole normalized I'}]
    return []

def load_current():
    manifest=_json(MANIFEST);_header(manifest)
    for section in ['files','source_bindings','code_bindings']:
        for rel,digest in manifest[section].items():
            if _sha(_path(ROOT,rel))!=digest:raise ValueError('Current '+section+' changed: '+rel)
    catalog=_bound(ROOT,manifest['catalog'])
    if catalog['physical_convention']!=CONVENTION or catalog['family_sign_conversion_allowed'] is not False:
        raise ValueError('Mixed current catalog convention')
    ids=[(x['loop_order'],x['orbit_id']) for x in catalog['entries']]
    if len(ids)!=len(set(ids)) or len(ids)!=163:raise ValueError('Current catalog orbit coverage mismatch')
    packages={}
    for ref in manifest['prepared_packages']:
        packages[ref['path']]=_bound(ROOT,ref);_native(packages[ref['path']])
    actual_counts={str(n):sum(x['loop_order']==n for x in catalog['entries']) for n in [1,2,3]}
    if actual_counts!={'1':2,'2':12,'3':149} or catalog.get('loop_counts')!=actual_counts:
        raise ValueError('Current catalog loop coverage mismatch')
    for row in catalog['entries']:
        if row['assembly_sign']!=1 or row['physical_convention']!=CONVENTION:raise ValueError('Current catalog sign mismatch')
        ref=row['native_prepared']
        if ref['path'] not in packages or ref['sha256']!=manifest['files'].get(ref['path']):raise ValueError('Unbound prepared orbit reference')
        _entry_definition(row,packages[ref['path']])
    return manifest,catalog

def check_current():
    manifest,catalog=load_current()
    return {'status':'passed','physical_convention':CONVENTION,'orbits':len(catalog['entries']),
      'loop_counts':catalog['loop_counts'],'prepared_packages':len(manifest['prepared_packages']),
      'bound_current_files':len(manifest['files']),'bound_immutable_sources':len(manifest['source_bindings']),
      'extra_family_signs':False,'formal':'not_run','scope':'Version/contract/coverage checks; use replay for actual assembly and RG checks.'}

def read_result(loop_order):
    manifest,_=load_current()
    return _bound(ROOT,manifest['results'][str(loop_order)])

def read_integral(orbit_id,definition=False,provenance=False):
    manifest,catalog=load_current();matches=[x for x in catalog['entries'] if x['orbit_id']==orbit_id]
    if len(matches)!=1:raise ValueError('Unknown or ambiguous current orbit ID')
    row=matches[0]
    if definition:
        data=_bound(ROOT,row['native_prepared']);value=_entry_definition(row,data)
        # Return the mathematical definition and current semantics without historical revision metadata.
        result={'orbit_id':orbit_id,'loop_order':row['loop_order'],'physical_convention':CONVENTION,
          'stage':'native full-color prepared C/J, not automatically a connected-I derivation',
          'C':value['C'],'J':value['J'],'normalization':value['normalization'],
          'preparation':value['preparation'],'quantity_boundary':row['quantity_boundary']}
    elif row.get('algebraic_zero'):
        result={'orbit_id':orbit_id,'loop_order':row['loop_order'],'physical_convention':CONVENTION,
                'status':'algebraic_zero','quantity_boundary':row['quantity_boundary'],'formal':'not_run'}
    elif row['loop_order']<=2:
        ref=row['current_coefficients'];result=_bound(ROOT,ref)
    else:
        selected=row['selected_evidence'];data=_bound(ROOT,selected['result'])
        # Distinct exact function/channel views remain distinct. No signs, values,
        # missing orders or representation identities are inferred by this reader.
        fields=_whole_i_views(selected['coefficient_field_candidates'],data)
        views=[{'quantity_field':x['field'],'quantity_role':x['quantity_role'],
                **({'channel':x['channel']} if 'channel' in x else {}),
                'powers':x['powers'],'coefficients':_pointer(data,x['field'])} for x in fields]
        result={'orbit_id':orbit_id,'loop_order':row['loop_order'],'physical_convention':CONVENTION,
          'normalization':selected['normalization'],'coefficient_views':views,
          'quantity_boundary':row['quantity_boundary'],'formal':'not_run',
          'status':'declared_whole_I_fields' if views else 'unavailable_use_joint_object_or_explicit_evidence'}
    if provenance:result=copy.deepcopy(result);result['provenance']=row
    return result

def replay_current(output_dir):
    manifest,_=load_current();out=Path(output_dir).resolve()
    if out.exists():raise ValueError('Use a new output directory')
    if out.is_relative_to(ROOT) and not any(out.is_relative_to(x) for x in [ROOT/'.wilson',ROOT/'replay-output']):
        raise ValueError('Inside the repository use .wilson or replay-output, never source/evidence trees')
    out.mkdir(parents=True);runs=[]
    for name in ['low','matter','AA','combine']:
        entry=manifest['replay'][name];dest=out/name
        cmd=[sys.executable,'-X','utf8','-B',str(_path(ROOT,entry['path'])),'--output-dir',str(dest)]
        if name=='matter':cmd+=['--lower-input',str(out/'low/coefficients.json')]
        if name=='combine':cmd+=['--low',str(out/'low/coefficients.json'),'--matter',str(out/'matter/matter.json'),'--aa',str(out/'AA/assembly.json')]
        p=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
        runs.append({'stage':name,'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
        (out/'runs.json').write_text(json.dumps(runs,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        if p.returncode:raise RuntimeError('Current '+name+' replay failed; see '+str(out/'runs.json'))
    load_current() # Refuse an input/source change during the replay.
    return {'status':'passed','output_dir':str(out),'stages':[x['stage'] for x in runs],
            'legacy_sign_conversion':False,'scope':'New direct coefficient assembly and RG; original integrations and dictionaries retain their stated reuse boundaries.','formal':'not_run'}

def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='command',required=True)
    sub.add_parser('check');p=sub.add_parser('result');p.add_argument('loop_order',type=int,choices=[1,2,3])
    p=sub.add_parser('integral');p.add_argument('orbit_id');p.add_argument('--definition',action='store_true');p.add_argument('--provenance',action='store_true')
    p=sub.add_parser('replay');p.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args()
    if a.command=='check':value=check_current()
    elif a.command=='result':value=read_result(a.loop_order)
    elif a.command=='integral':value=read_integral(a.orbit_id,a.definition,a.provenance)
    else:value=replay_current(a.output_dir)
    print(json.dumps(value,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
