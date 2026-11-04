"""Exact new printed-function/operator transcription checks; no integration."""
from pathlib import Path
import importlib.util,json,hashlib,re
import sympy as s
from sympy.parsing.sympy_parser import parse_expr

ROOT=Path(__file__).resolve().parents[4];HERE=Path(__file__).resolve().parent
FAMILY=ROOT/'integrals/families/three-loop-insertions';SNAP=FAMILY/'runtime/snapshots'
PARSER=ROOT/'workflow/checks/integral-documentation-6e-20261103/review/check_hlog_transcription.py'
spec=importlib.util.spec_from_file_location('reviewed_hlog_transcription_parser',PARSER);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def text_expression(value):
    if '=' in value:value=value.split('=',1)[1]
    value=re.sub(r'Hlog\(z,\s*[\[(]([^\])]*?)[\])]\)',lambda m:h.symbol(tuple(int(v) for v in m[1].split(',') if v.strip())),value)
    value=re.sub(r'zeta\[(\d+)\]',lambda m:'Z'+m[1],value).replace('^','**')
    return s.expand(parse_expr(value,local_dict={**h.LOCALS,'zeta':s.zeta}))

def tex_expression(value):
    h.LOCALS['e']=s.Symbol('e')
    value=re.sub(r'(?<=\d)e',r'*e',value) # avoid scientific-number lexing of 2e-3
    def word(m):
        raw=m[1] if m[1] is not None else m[2]
        letters=tuple(map(int,raw.split(','))) if ',' in raw or raw=='-1' else (0,)*len(raw)
        return '(G('+','.join(map(str,letters))+';z))'
    value=re.sub(r'G_\{([^}]+)\}|G_(-?\d+)',word,value)
    value=re.sub(r'\\zeta_(\d+)',lambda m:r'\zeta_{'+m[1]+'}',value)
    value=value.replace('\\\\','').replace('&','').replace(r'\tfrac',r'\frac').replace('{}','')
    value=re.sub(r'\\frac(\{[^{}]*\})([0-9])',lambda m:r'\frac'+m[1]+'{'+m[2]+'}',value)
    value=re.sub(r'\}\s+\{','}{',value)
    return h.latex(value.replace('\n',' '))

def record(label,md,old_path,a,b):
    delta=s.cancel(s.expand(a-b));assert delta==0,(label,delta)
    symbols={v for v in a.free_symbols if str(v) in h.WORDS}
    return {'label':label,'status':'passed','new_markdown':md.relative_to(ROOT).as_posix(),'new_markdown_sha256':sha(md),'source':old_path.relative_to(ROOT).as_posix(),'source_sha256':sha(old_path),'nonzero_HLOG_coordinates':len(symbols)+(a.subs({v:0 for v in symbols})!=0)}

def main():
    rows=[];md=FAMILY/'appendices/pi1-moments.md';body=md.read_text(encoding='utf-8')
    blocks=re.findall(r'```text\s*\n(.*?)\n```',body,re.S);assert len(blocks)==11,len(blocks)
    expected={}
    for block in blocks:
        name,value=block.split('=',1);name=name.strip()
        if name.startswith('M_'):
            path=SNAP/'workflow/checks/single-opposite/runs/20260916T090538777022Z'/(name+'.txt')
            old=text_expression(path.read_text(encoding='utf-8').strip().rstrip(':;'))
        else:
            assert name=='rho1';path=SNAP/'benchmarks/2ne/families/same-channel-two/structural-basis-functions.json'
            old=text_expression(json.loads(path.read_text())['functions']['rho_id_0_0_1']['hlog_expression'])
        actual=text_expression(value);rows.append(record(name,md,path,actual,old));expected[name]=old
    md=FAMILY/'chains-pi2.md';body=md.read_text(encoding='utf-8')
    patterns={'M_0_1':r'M_\{0,1\}=([^\n]+?),\\qquad',
              'M_1_1':r'M_\{1,1\}=([^\n]+),\s*\n\\\]',
              'M_2_1':r'M_\{2,1\}=\{\}&(.*?),\\\\\s*\nM_\{0,3\}',
              'M_0_3':r'M_\{0,3\}=\{\}&(.*?)=-M_\{0,1\}\^2'}
    for name,pattern in patterns.items():
        match=re.search(pattern,body,re.S);assert match,(name,'not found')
        path=SNAP/'workflow/checks/single-opposite/runs/20260916T090538777022Z'/(name+'.txt')
        rows.append(record(name+' printed LaTeX',md,path,tex_expression(match[1]),expected[name]))
    path=SNAP/'benchmarks/current/prepared/3loop/2e-Pi2.json';native=json.loads(path.read_text())['orbits'][0]['operator']['channels']
    for key in ('A','B'):
        for channel in ('AA','Af','Ff'):
            pattern=re.escape(key+'_{'+channel+'}')+r'&=(.*?)(?:,&|,\\\\|\.\s*\n\\end)'
            match=re.search(pattern,body,re.S);assert match,(key,channel)
            value=match[1]
            if channel=='Ff' and key=='B':value=value.rstrip('.')
            e=s.Symbol('e');old=s.sympify(native[channel][key].replace('^','**'),locals={'eps':e})
            rows.append(record(key+'_'+channel,md,path,tex_expression(value),old))
    result={'status':'passed','expressions':len(rows),'HLOG_nonzero_coordinate_positions':sum(int(r['nonzero_HLOG_coordinates']) for r in rows[:15]),'parser_source':PARSER.relative_to(ROOT).as_posix(),'parser_sha256':sha(PARSER),'checks':rows,'scope':'Actual new Markdown native/LaTeX expressions versus exact frozen adopted M/rho1 evaluations and native Pi2 A/B rational fields. This verifies mathematical transcription, not the original integrations or the literature operator.'}
    (HERE/'math-transcription-check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print({k:result[k] for k in ('status','expressions','HLOG_nonzero_coordinate_positions')})
if __name__=='__main__':main()
