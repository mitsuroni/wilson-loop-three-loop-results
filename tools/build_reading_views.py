"""Use GitHub math delimiters in reader pages; preserve frozen result bytes.

Run from any directory with Python 3.12+. No network or scientific CAS is used.
"""
from pathlib import Path
import hashlib, json, os, re
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
PACKAGE=ROOT/'results/2026-10-23'
GUIDE=ROOT/'docs/data-flow'
VIEWS=ROOT/'docs/results'
CODE=re.compile(r'\$`[^\n]*?`\$|```[\s\S]*?```|`[^`\n]+`')
MATH=re.compile(r'\\\[([\s\S]*?)\\\]|\\\(([^\n]*?)\\\)')
ANY_MATH=re.compile(r'```math\n([\s\S]*?)\n```|\\\[([\s\S]*?)\\\]|\\\(([^\n]*?)\\\)|\$\$([\s\S]*?)\$\$|\$`([^\n]*?)`\$|(?<!\\)\$(?!\$)([^\n]*?)(?<!\\)\$')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def math_bodies(text):
    text=CODE.sub(lambda m:m.group() if m.group().startswith(('```math\n','$`')) else '',text)
    return [(('block' if m.group(1) is not None or m.group(2) is not None or m.group(4) is not None else 'inline'),next(x for x in m.groups() if x is not None)) for m in ANY_MATH.finditer(text)]

def convert(text):
    """Change only math delimiters outside code; preserve each TeX body."""
    saved=[]
    def protect(m):
        if m.group().startswith(('```math\n','$`')):return m.group()
        saved.append(m.group());return f'CODELITERAL{len(saved)-1}END'
    protected=CODE.sub(protect,text)
    changes=[]
    def math(m):
        block=m.group(1) is not None or m.group(2) is not None or m.group(4) is not None
        body=next(x for x in m.groups() if x is not None)
        if m.group(2) is not None or m.group(3) is not None:
            changes.append({'kind':'block' if block else 'inline','tex_sha256':hashlib.sha256(body.encode()).hexdigest()})
        if not block:return '$`'+body+'`$'
        # GitHub math fences preserve TeX backslashes through Markdown parsing.
        before,after=protected[:m.start()],protected[m.end():]
        prefix='' if not before or before.endswith('\n\n') else ('\n' if before.endswith('\n') else '\n\n')
        suffix='' if not after or after.startswith('\n\n') else ('\n' if after.startswith('\n') else '\n\n')
        return prefix+'```math\n'+body+'\n```'+suffix
    result=ANY_MATH.sub(math,protected)
    assert not MATH.search(result)
    for i,code in enumerate(saved):result=result.replace(f'CODELITERAL{i}END',code)
    assert math_bodies(text)==math_bodies(result),'A mathematical expression changed'
    return result,changes

def rewrite_links(text,src,dst):
    def target(href):
        if href.startswith(('http:','https:','#','mailto:')):return href
        raw,sep,frag=href.partition('#')
        resolved=(src.parent/unquote(raw)).resolve()
        if resolved.parent==PACKAGE and resolved.suffix=='.md':resolved=VIEWS/resolved.name
        return os.path.relpath(resolved,dst.parent).replace('\\','/')+(sep+frag if sep else '')
    # Only prose links; source formulas inside fenced/inline code stay literal.
    pieces=CODE.split(text);spans=list(CODE.finditer(text));out=[]
    for i,piece in enumerate(pieces):
        piece=re.sub(r'\]\(([^)\n]+)\)',lambda m:']('+target(m.group(1))+')',piece)
        piece=re.sub(r'(src=")([^"]+)(")',lambda m:m.group(1)+target(m.group(2))+m.group(3),piece)
        out.append(piece)
        if i<len(spans):out.append(spans[i].group())
    return ''.join(out)

def main():
    VIEWS.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((PACKAGE/'package-manifest.json').read_text(encoding='utf-8'))['files']
    for rel,info in manifest.items():assert sha(PACKAGE/rel)==info['sha256'],rel
    result={'format':'wilson-reader-math-repair/1','date':'2026-10-23','preserved_package_files':len(manifest),'pages':{}}
    audit_path=ROOT/'docs/reading-view-manifest.json'
    previous=json.loads(audit_path.read_text(encoding='utf-8')) if audit_path.exists() else {'pages':{}}
    for src in sorted(PACKAGE.glob('*.md')):
        dst=VIEWS/src.name
        text,changes=convert(src.read_text(encoding='utf-8'))
        text=rewrite_links(text,src,dst)
        dst.write_text(text,encoding='utf-8',newline='\n')
        result['pages'][dst.relative_to(ROOT).as_posix()]={'canonical_source':src.relative_to(ROOT).as_posix(),'source_sha256':sha(src),'converted_math':changes}
    for src in [GUIDE/'README.md',GUIDE/'reading-guide-zh.md']:
        old=src.read_text(encoding='utf-8');text,changes=convert(old)
        text=rewrite_links(text,src,src)
        src.write_text(text,encoding='utf-8',newline='\n')
        # Refresh only the embedded prose; keep the existing offline MathJax assets.
        htmlfile=GUIDE/('index.html' if src.name=='README.md' else 'reading-guide-zh.html')
        rendered=text
        for name in ['README','reading-guide-zh','source-register']:
            rendered=rendered.replace(']('+name+'.md)',']('+('index' if name=='README' else name)+'.html)')
        # Result reading views are for GitHub; the offline links use frozen HTML.
        rendered=re.sub(r'\]\(../results/([^)]+)\.md\)',r'](../../results/2026-10-23/\1.html)',rendered)
        rendered=re.sub(r'```math\n([\s\S]*?)\n```',lambda m:'\\['+m.group(1)+'\\]',rendered)
        rendered=re.sub(r'\$`([^\n]*?)`\$',lambda m:'\\('+m.group(1)+'\\)',rendered)
        data=json.dumps(rendered,ensure_ascii=False).replace('</',r'<\/')
        oldhtml=htmlfile.read_text(encoding='utf-8')
        newhtml,count=re.subn(r'const markdown=[\s\S]*?;\s*\nconst saved=',lambda m:'const markdown='+data+';\nconst saved=',oldhtml,count=1)
        assert count==1
        htmlfile.write_text(newhtml,encoding='utf-8',newline='\n')
        rel=src.relative_to(ROOT).as_posix()
        if not changes:changes=previous['pages'].get(rel,{}).get('converted_math',[])
        result['pages'][rel]={'converted_math':changes,'html_refreshed':True}
    readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8')
    text=text.replace('(results/2026-10-23/index.md)','(docs/results/index.md)').replace('(results/2026-10-23/conventions.md)','(docs/results/conventions.md)').replace('(results/2026-10-23/validation-scope.md)','(docs/results/validation-scope.md)')
    readme.write_text(text,encoding='utf-8',newline='\n')
    links=0
    for src in [GUIDE/'README.md',GUIDE/'reading-guide-zh.md',*VIEWS.glob('*.md')]:
        prose=CODE.sub('',src.read_text(encoding='utf-8'))
        assert not MATH.search(prose),src
        for href in re.findall(r'\]\(([^)\n]+)\)',prose):
            if href.startswith(('http:','https:','#')):continue
            assert (src.parent/unquote(href.split('#')[0])).exists(),(src,href)
            links+=1
    for rel,info in manifest.items():assert sha(PACKAGE/rel)==info['sha256'],rel
    result['local_links']=links
    result['scientific_data_unchanged']=True
    # Audit the output as rendered source, rather than writing volatile run paths.
    (ROOT/'docs/reading-view-manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'pages':len(result['pages']),'links':links,'preserved_files':len(manifest),'conversions':sum(len(p['converted_math']) for p in result['pages'].values())}))

if __name__=='__main__':main()
