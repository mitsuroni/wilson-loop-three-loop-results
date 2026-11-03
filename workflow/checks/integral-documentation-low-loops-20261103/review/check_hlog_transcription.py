"""Check new mathematical tables against adopted raw jets; no integration."""
from pathlib import Path
import hashlib
import json
import re
import sympy as s
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
MD=ROOT/'integrals/two-loop/3e-3g-04/appendices/hlog-jets.md'
SOURCE=ROOT/'workflow/checks/tgv04-high-order/runs/positive-epsilon'
TRANSFORMS=standard_transformations+(implicit_multiplication_application,)
LOCALS={'z':s.Symbol('z'),'Z2':s.Symbol('Z2'),'Z3':s.Symbol('Z3'),'Z5':s.Symbol('Z5')}

def group(text,pos):
    assert text[pos]=='{',text[pos:]
    depth=1
    end=pos+1
    while depth:
        if text[end]=='{':depth+=1
        if text[end]=='}':depth-=1
        end+=1
    return text[pos+1:end-1],end

def latex(text):
    text=text.strip().strip('$').replace(r'\left','').replace(r'\right','')
    text=re.sub(r'\\zeta_\{(\d+)\}',lambda m:'Z'+m[1],text)
    while r'\frac' in text:
        start=text.rfind(r'\frac')
        a,end=group(text,start+5)
        b,end=group(text,end)
        text=text[:start]+'(('+a+')/('+b+'))'+text[end:]
    text=text.replace('^','**').replace('{','(').replace('}',')')
    return parse_expr(text,local_dict=LOCALS,transformations=TRANSFORMS)

def native(path):
    text=path.read_text(encoding='utf-8').split('=',1)[1].strip().rstrip(':;')
    symbols={}
    def hlog(m):
        word=tuple(int(v.strip()) for v in m[1].split(',') if v.strip())
        name='HH'+str(len(symbols))
        symbols[word]=s.Symbol(name)
        return name
    text=re.sub(r'Hlog\(z,\[([^\]]*)\]\)',hlog,text)
    text=re.sub(r'zeta\[(\d+)\]',lambda m:'Z'+m[1],text).replace('^','**')
    local=dict(LOCALS);local.update({str(v):v for v in symbols.values()})
    expression=s.expand(parse_expr(text,local_dict=local))
    output={w:expression.coeff(v) for w,v in symbols.items()}
    output[()]=expression.subs({v:0 for v in symbols.values()})
    return {w:c for w,c in output.items() if c!=0}

def main():
    tables={};current=None
    for line in MD.read_text(encoding='utf-8').splitlines():
        section=re.fullmatch(r'## \$([DHO])_(\d)\$',line)
        if section:current=section[1]+section[2];tables[current]={};continue
        if current and line.startswith('## '):current=None
        if not current or not line.startswith('| $'):continue
        cells=[v.strip().strip('$') for v in line.strip('|').split('|')]
        word=() if cells[0]==r'\varnothing' else tuple(int(v) for v in cells[0].strip('()').split(','))
        assert word not in tables[current],(current,word)
        tables[current][word]=latex(cells[1])
    checks=[]
    for name,table in tables.items():
        source=SOURCE/(name+'.txt');old=native(source)
        for word in set(old)|set(table):
            delta=s.cancel(old.get(word,0)-table.get(word,0))
            assert delta==0,(name,word,delta)
        checks.append({'jet':name,'all_coordinates_equal':True,'coordinates':len(table),'source':str(source.relative_to(ROOT)).replace('\\','/'),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
    assert set(tables)=={'D2','D3','H1','H2','O1','O2'}
    result={'status':'passed','scope':'Exact new Markdown HLOG-coordinate transcription versus adopted original raw jet outputs. No integration, no independent validation of the original HyperInt evaluations.','markdown':str(MD.relative_to(ROOT)).replace('\\','/'),'markdown_sha256':hashlib.sha256(MD.read_bytes()).hexdigest(),'jets':checks}
    (HERE/'hlog-transcription-check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','jets':len(checks),'exact_coordinate_checks':sum(v['coordinates'] for v in checks)}))

if __name__=='__main__':main()
