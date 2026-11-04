"""Render review generators without changing historical review/evidence files."""
import ast,contextlib,io,json,re,runpy
from pathlib import Path
from unittest.mock import patch
from update import ROOT,HERE,BATCHES,LABELS,sha,put
def main():
    results=[]
    for name in BATCHES:
        base=ROOT/'workflow/checks'/name/'review';script=base/'write_review.py';ast.parse(script.read_text(encoding='utf-8-sig'));captured={}
        def capture(p,data,*args,**kwargs):
            p=p.resolve();assert p in (base/'audit.json',base/'README.md'),p;captured[p]=data;return len(data)
        with patch.object(Path,'write_text',capture),contextlib.redirect_stdout(io.StringIO()):runpy.run_path(str(script),run_name='__main__')
        # The historical writer can observe subsequently updated source hashes.
        # Never replace the scientific audit with that fresh reconstruction.
        # Compare the actual edited record to its preserved historical bytes.
        previous=json.loads((HERE/'before'/name/'audit.json').read_text(encoding='utf-8-sig'))
        current=json.loads((base/'audit.json').read_text(encoding='utf-8-sig'))
        for row in previous['entries']:row['document_status']=LABELS[row['document_status']]
        assert previous==current,'Actual mathematical review fields changed: '+name
        assert not re.search('[\u4e00-\u9fff]',captured[base/'README.md']),name
        actual=(base/'README.md').read_text(encoding='utf-8-sig')
        assert not re.search('[\u4e00-\u9fff]',actual),name
        # Current and regenerated tables must have identical labels and rows.
        assert [s for s in actual.splitlines() if s.startswith('|')]==[s for s in captured[base/'README.md'].splitlines() if s.startswith('|')],name
        for p,data in captured.items():
            dst=HERE/'generated'/name/p.name;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_text(data,encoding='utf-8')
        results.append({'batch':name,'entries':len(json.loads(captured[base/'audit.json'])['entries']),'regenerated_English_only':True,'audit_json_unchanged_except_prior_label_translation':True,'rendered_table_matches_generator':True})
    # Check all current peer-facing Markdown in these three review directories.
    pages=[];links=0
    for name in BATCHES:
        for p in (ROOT/'workflow/checks'/name/'review').glob('*.md'):
            text=p.read_text(encoding='utf-8-sig');assert not re.search('[\u4e00-\u9fff]',text),p
            for m in re.finditer(r'\[[^\]\n]+\]\(([^)\n]+)\)',text):
                url=m.group(1)
                if '://' in url or url.startswith('#'):continue
                assert (p.parent/url.split('#')[0]).resolve().exists(),(p,url);links+=1
            pages.append(p.relative_to(ROOT).as_posix())
    put(HERE/'check.json',{'status':'passed','review_entries':sum(r['entries'] for r in results),'review_pages_checked':len(pages),'local_links_checked':links,'generators':results,'scope':'Language, generator/table agreement and exact audit-field preservation only. No mathematical evaluation or new review decision.'})
    print({'status':'passed','entries':68,'English_review_pages':len(pages),'generator_regressions':0})
if __name__=='__main__':main()
