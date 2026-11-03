"""Normalize trusted inputs before Wilson; no Lean, MCP, receipts, or database access."""
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
from .color import project,A,F,require
from .ordered import REPRESENTATION, integral_ast

ROOT=Path(__file__).resolve().parents[1]
CORPUS=ROOT/'benchmarks/2ne'
PREPARATION_VERSION='ordered-color-euler-input/0.3'
MANIFEST_FORMAT='wilson-prepared-input-manifest/0.3'
RECIPE='ordered-parameters+fixed-loop-color-a1+per-loop-euler-compensation/0.3'
EULER_SCHEME='exp(-loop_order*EulerGamma*eps)'

def encoded(data):return (json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+'\n').encode('utf-8')
def digest(raw):return hashlib.sha256(raw).hexdigest()
def display(expr):return str(expr).replace('**','^')

def euler_factor(n):
    return 'Exp[-EulerGamma*eps]' if n==1 else f'Exp[-{n}*EulerGamma*eps]'

def validate_package_normalization_input(package):
    require(package.get('package_kind')=='montopar-initial-integrals','unexpected source format')
    n=package.get('loop_order')
    require(isinstance(n,int) and not isinstance(n,bool) and n>0,'loop_order must be a positive integer')
    require(package.get('emission_type')==f'{2*n}e','emission family does not match loop_order')
    require(isinstance(package.get('orbits'),list),'missing orbit list')
    expected_ft=f'Gamma[1-eps]^{n}/(4*Pi^(2-eps))^{n}'
    for orbit in package['orbits']:
        integral_ast(orbit['J'])
        require(isinstance(orbit.get('C'),dict),'missing C')
        combined=orbit['C'].get('combined')
        components=orbit['C'].get('components')
        require(isinstance(combined,str) and isinstance(components,dict),'missing C expression or components')
        require('euler_compensation' not in components,'Euler compensation is already declared')
        require('EulerGamma' not in combined and all('EulerGamma' not in str(v) for v in components.values()),
            'Euler compensation is already present')
        require(components.get('ft_prefactor')==expected_ft,'Gamma power does not match loop_order')
        require(f'({expected_ft})' in combined,'C.combined does not contain its declared Fourier/Gamma factor')
    return n

def atomic_write(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and path.read_bytes()==raw:return
    fd,name=tempfile.mkstemp(dir=path.parent,prefix='.prepare-')
    try:
        with os.fdopen(fd,'wb') as stream:stream.write(raw)
        os.replace(name,path)
    finally:
        Path(name).unlink(missing_ok=True)

def normalize_package(package,source,cache):
    n=validate_package_normalization_input(package)
    output=copy.deepcopy(package)
    output['package_kind']='wilson-prepared-integrals'
    output['preparation']=dict(version=PREPARATION_VERSION,source=source,
        source_package_kind=package['package_kind'],trust_basis='user_accepted_preprocessed_input',
        proof_boundary='Wilson verification begins with prepared C and the pre-hypercube J with its explicit iterated order.',
        integral_representation=REPRESENTATION,
        color_policy='fixed-loop full coordinates, select a1',
        other_C_factors='unchanged except the declared Euler compensation factor',
        euler_compensation=EULER_SCHEME,lean_invoked=False)
    compensation=euler_factor(n)
    for orbit in output['orbits']:
        original=copy.deepcopy(orbit['C']);text=original['components']['color'];key=(n,text)
        if key not in cache:cache[key]=project(n,text)
        result=cache[key];coordinate=result['coordinates'][0]
        selected=display(coordinate*A**(n-1)*F)
        prefix='('+text+')'
        require(original['combined'].startswith(prefix+' * '),'C.combined does not begin with its declared color factor')
        orbit['C']['components']['color']=selected
        orbit['C']['components']['euler_compensation']=compensation
        orbit['C']['combined']='('+selected+') * ('+compensation+')'+original['combined'][len(prefix):]
        orbit['color_preparation']=dict(original_C=original,
            basis=[display(b) for b in result['basis']],
            coordinates=[str(c) for c in result['coordinates']],a1=str(coordinate),
            reconstructed_color=display(result['reconstruction']),projected_color=selected,
            python_checks=dict(unique_rational_coordinates=True,exact_reconstruction_residual_zero=True),
            evidence_kind='python_exact_preprocessing',formal_receipt=None)
        orbit['euler_preparation']=dict(loop_order=n,gamma_power=n,
            original_ft_prefactor=original['components']['ft_prefactor'],factor=compensation,
            scheme=EULER_SCHEME,
            python_checks=dict(loop_order_positive=True,emission_family_matches=True,
                gamma_power_matches_loop_order=True,previous_compensation_absent=True),
            evidence_kind='python_deterministic_preprocessing',formal_receipt=None)
    return output

def verify_prepared(output=CORPUS/'prepared'):
    output=Path(output).resolve()
    require(output.is_relative_to(ROOT),'prepared directory must remain in the workspace')
    manifest=json.loads((output/'manifest.json').read_text(encoding='utf-8'))
    require(manifest.get('format')==MANIFEST_FORMAT,'unexpected prepared manifest format')
    require(manifest.get('recipe')==RECIPE,'unexpected prepared recipe')
    require(manifest.get('euler_compensation')==EULER_SCHEME,'unexpected Euler compensation scheme')
    expected_code={p.relative_to(ROOT).as_posix():digest(p.read_bytes())
        for p in [Path(__file__),Path(__file__).with_name('color.py'),Path(__file__).with_name('ordered.py')]}
    require(manifest.get('code_sha256')==expected_code,'preprocessing code checksum mismatch')
    total=0
    for entry in manifest.get('packages',[]):
        path=(output/entry['path']).resolve()
        require(path.parent==output,'prepared package path outside output directory')
        raw=path.read_bytes()
        require(digest(raw)==entry['sha256'],'prepared package checksum mismatch')
        package=json.loads(raw.decode('utf-8'))
        require(package.get('package_kind')=='wilson-prepared-integrals','unexpected prepared package format')
        require(package.get('preparation',{}).get('version')==PREPARATION_VERSION,'prepared package version mismatch')
        require(package.get('preparation',{}).get('euler_compensation')==EULER_SCHEME,'prepared package compensation mismatch')
        require(package.get('preparation',{}).get('integral_representation')==REPRESENTATION,'prepared representation mismatch')
        for orbit in package['orbits']: integral_ast(orbit['J'])
        require(package.get('emission_type')==entry['family'],'prepared family mismatch')
        ids=[orbit.get('orbit_id') for orbit in package.get('orbits',[])]
        require(ids==entry['orbit_ids'] and len(ids)==entry['count'],'prepared orbit inventory mismatch')
        total+=entry['count']
    require(total==manifest.get('total_orbits'),'prepared total orbit count mismatch')
    return manifest

def prepare(corpus=CORPUS,output=None):
    corpus=Path(corpus).resolve();output=Path(output or corpus/'prepared').resolve()
    require(output.is_relative_to(ROOT),'output must remain in the workspace')
    require(not output.is_relative_to(corpus/'sources'),'never overwrite trusted raw sources')
    require(output!=corpus,'use a separate output directory')
    manifest=json.loads((corpus/'ordered-source-manifest.json').read_text(encoding='utf-8'))
    require(manifest.get('representation')==REPRESENTATION,'ordered source representation required')
    cache={};pending=[];entries=[]
    for item in manifest['sources']:
        path=(corpus/item['path']).resolve()
        require(path.is_relative_to(corpus),'source path outside corpus')
        raw=path.read_bytes();require(digest(raw)==item['sha256'],'raw source checksum mismatch')
        package=json.loads(raw.decode('utf-8-sig'))
        require(package['package_kind']=='montopar-initial-integrals','unexpected source format')
        require(package['emission_type']==item['family'] and len(package['orbits'])==item['orbit_count'],'source inventory mismatch')
        ids=[o['orbit_id'] for o in package['orbits']]
        require(ids==item['orbit_ids'] and len(ids)==len(set(ids)),'source orbit IDs mismatch')
        source=dict(path_base=corpus.relative_to(ROOT).as_posix(),path=item['path'],sha256=item['sha256'])
        prepared=normalize_package(package,source,cache);content=encoded(prepared)
        filename=item['family']+'.json'
        require(Path(filename).name==filename,'invalid output filename')
        pending.append((output/filename,content))
        entries.append(dict(family=item['family'],path=filename,sha256=digest(content),source=source,orbit_ids=ids,count=len(ids)))
    index=dict(format=MANIFEST_FORMAT,trust_basis='user_accepted_preprocessed_input',
        recipe=RECIPE,euler_compensation=EULER_SCHEME,integral_representation=REPRESENTATION,lean_invoked=False,
        total_orbits=sum(x['count'] for x in entries),unique_color_expressions=len(cache),packages=entries,
        code_sha256={p.relative_to(ROOT).as_posix():digest(p.read_bytes()) for p in [Path(__file__),Path(__file__).with_name('color.py'),Path(__file__).with_name('ordered.py')]})
    # Compute/validate every package before publishing; manifest is published last.
    for path,raw in pending:atomic_write(path,raw)
    atomic_write(output/'manifest.json',encoded(index))
    return verify_prepared(output)
