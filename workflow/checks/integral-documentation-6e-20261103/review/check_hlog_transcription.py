"""Exact new HLOG prose/table transcription checks. No integrals are evaluated."""
from pathlib import Path
import hashlib
import json
import re
import sympy as s
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
TRANSFORMS=standard_transformations+(implicit_multiplication_application,)
LOCALS={'z':s.Symbol('z'),'pi':s.pi}
LOCALS.update({'Z'+str(n):s.zeta(n) for n in range(2,9)})
WORDS={}

def symbol(word):
    if not word:return '1'
    name='HH'+'_'.join('m1' if a==-1 else str(a) for a in word)
    WORDS[name]=tuple(word);LOCALS[name]=s.Symbol(name)
    return name

def group(text,pos):
    assert text[pos]=='{',text[pos:]
    depth=1;end=pos+1
    while depth:
        if text[end]=='{':depth+=1
        if text[end]=='}':depth-=1
        end+=1
    return text[pos+1:end-1],end

def latex(text):
    text=text.strip().strip('$').rstrip('.').replace(r'\left','').replace(r'\right','')
    text=re.sub(r'G\(([-0-9,]*);z\)',lambda m:symbol(tuple(int(v) for v in m[1].split(',') if v)),text)
    text=re.sub(r'\\zeta\^\{(\d+)\}\((\d+)\)',lambda m:'(Z'+m[2]+'^'+m[1]+')',text)
    text=re.sub(r'\\zeta_\{(\d+)\}',lambda m:'Z'+m[1],text)
    text=re.sub(r'\\zeta\((\d+)\)',lambda m:'Z'+m[1],text)
    while r'\frac' in text:
        start=text.rfind(r'\frac');a,end=group(text,start+5);b,end=group(text,end)
        text=text[:start]+'(('+a+')/('+b+'))'+text[end:]
    text=text.replace(r'\pi','pi').replace(r'\cdot','*').replace('^','**').replace('{','(').replace('}',')')
    if '\\' in text:raise ValueError('Unparsed LaTeX: '+text)
    return parse_expr(text,local_dict=LOCALS,transformations=TRANSFORMS)

def native(path):
    text=path.read_text(encoding='utf-8').strip().rstrip(':;')
    if re.match(r'^\w+\s*=',text):text=text.split('=',1)[1]
    text=re.sub(r'Hlog\(z,\s*\[([^\]]*)\]\)',lambda m:symbol(tuple(int(v.strip()) for v in m[1].split(',') if v.strip())),text)
    text=re.sub(r'zeta\[(\d+)\]',lambda m:'Z'+m[1],text).replace('^','**')
    return s.expand(parse_expr(text,local_dict=LOCALS))

def bind(md,link):
    path=(md.parent/link).resolve()
    assert path.is_relative_to(ROOT/'integrals/families/three-loop-6e/runtime/primitives'),path
    assert path.is_file(),path
    return path

def compare(md,label,source,actual):
    old=native(source)
    delta=s.cancel(s.expand(old-actual))
    assert delta==0,(md,label,delta)
    symbols=sorted(old.free_symbols|actual.free_symbols,key=str)
    count=sum(1 for v in symbols if str(v) in WORDS)
    if old.subs({v:0 for v in symbols if str(v) in WORDS})!=0:count+=1
    return {'label':label,'status':'passed','exact_nonzero_coordinates':count,'source':str(source.relative_to(ROOT)).replace('\\','/'),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}

def check_single(md):
    source=None;records=[]
    for line in md.read_text(encoding='utf-8').splitlines():
        if line.startswith('Order'):
            source=bind(md,re.search(r'\]\(([^)]+)\)',line)[1])
        if line.startswith('$$h_'):
            label,expression=line.strip('$').split('=',1)
            records.append(compare(md,label,source,latex(expression.rstrip('$'))))
    assert len(records)==30,len(records)
    return records

def check_tables(md):
    label=None;source=None;value=s.S.Zero;records=[];have=False
    def finish():
        if label is not None:
            assert source is not None,(md,label)
            assert have,(md,label)
            records.append(compare(md,label,source,value))
    for line in md.read_text(encoding='utf-8').splitlines():
        if line.startswith('## $'):
            finish();label=line[3:].strip('$');source=None;value=s.S.Zero;have=False
        elif label and line.startswith('Source:'):
            source=bind(md,re.search(r'\]\(([^)]+)\)',line)[1])
        elif label and line.startswith('| $'):
            cells=[v.strip().strip('$') for v in line.strip('|').split('|')]
            word=() if cells[0]==r'\varnothing' else tuple(int(v) for v in cells[0].strip('()').split(','))
            value+=latex(cells[1])*s.sympify(symbol(word),locals=LOCALS);have=True
        elif label and line.startswith('$$') and re.search(r'=0\.\$\$$',line):
            have=True
    finish()
    expected=2 if 'original zero-order' in md.read_text(encoding='utf-8').splitlines()[0] else 6
    assert len(records)==expected,(md,len(records),expected)
    return records

def main():
    files=[ROOT/'integrals/families/three-loop-6e/appendices/single-opposite-jets.md']
    ids=[10,11,12,20,21,22,23,24,25,32,33,34,37,38]
    files += [ROOT/f'integrals/three-loop/6e-{i:02}/appendices/evaluated-jets.md' for i in ids]
    output=[]
    for md in files:
        assert md.is_file(),md
        records=check_single(md) if 'single-opposite-jets' in md.name else check_tables(md)
        output.append({'markdown':str(md.relative_to(ROOT)).replace('\\','/'),'markdown_sha256':hashlib.sha256(md.read_bytes()).hexdigest(),'records':records})
    result={'status':'passed','scope':'Exact mathematical transcription of new HLOG expressions/tables versus their declared immutable raw CAS outputs. No integration or independent verification of the original CAS evaluations.','files':output,'evaluated_expressions':sum(len(v['records']) for v in output),'exact_nonzero_coordinates':sum(r['exact_nonzero_coordinates'] for v in output for r in v['records'])}
    (HERE/'hlog-transcription-check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','evaluated_expressions','exact_nonzero_coordinates']}))

if __name__=='__main__':main()
