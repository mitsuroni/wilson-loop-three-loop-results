"""English-only presentation of existing reviews, preserving mathematical evidence."""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
BATCHES=('integral-documentation-low-loops-20261103','integral-documentation-6e-20261103','integral-documentation-insertions-20261104')
LABELS={'正文可连续阅读（接受明确的物理算符根）':'Complete','正文可连续阅读':'Complete','明确采用已求函数引理，链已接通':'Complete (adopted lemmas)','明确采用已求引理，链已接通':'Complete (adopted lemmas)'}
LEGEND='Status describes the completeness of the written derivation. "Adopted lemmas" means that explicitly stated prior evaluations are used; accepted physical inputs and proof limitations remain as documented.'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
    assert not (HERE/'changes.json').exists(),'This editorial revision has already run.'
    changes=[]
    for name in BATCHES:
        base=ROOT/'workflow/checks'/name/'review'
        for filename in ('README.md','write_review.py','audit.json'):
            p=base/filename;old=p.read_text(encoding='utf-8-sig');new=old
            backup=HERE/'before'/name/filename;backup.parent.mkdir(parents=True,exist_ok=True);backup.write_bytes(p.read_bytes())
            for a,b in LABELS.items():new=new.replace(a,b)
            if filename!='audit.json':
                new=re.sub(r'\b([Aa]ll|not)(?=(?:14|15|39)\b)',r'\1 ',new)
                if filename=='README.md':
                    at=new.index('\n| Entry |');new=new[:at]+'\n'+LEGEND+'\n'+new[at:]
                else:
                    # Place the same explanation in the generated Markdown.
                    header=next(s for s in re.findall(r"'([^'\n]*\| Entry \|[^'\n]*)'",new))
                    literal=repr(header)
                    assert literal in new,(p,header)
                    new=new.replace(literal,repr(LEGEND)+", '', "+literal,1)
            else:
                a=json.loads(old);b=json.loads(new)
                for previous,current in zip(a['entries'],b['entries']):
                    assert current['document_status']==LABELS[previous['document_status']]
                    current['document_status']=previous['document_status']
                assert a==b,'Only presentation labels may change in the audit.'
            p.write_text(new,encoding='utf-8');changes.append({'path':p.relative_to(ROOT).as_posix(),'before_sha256':sha(backup),'after_sha256':sha(p),'historical_copy':backup.relative_to(ROOT).as_posix()})
    # Translate the one quoted Chinese phrase in an otherwise English guide.
    p=ROOT/'workflow/checks'/BATCHES[0]/'review/standard-forward-test.md'
    old=p.read_text(encoding='utf-8-sig');new=old.replace('“按需” rule','requirement to include material only as needed')
    if new!=old:
        backup=HERE/'before'/BATCHES[0]/p.name;backup.write_bytes(p.read_bytes());p.write_text(new,encoding='utf-8');changes.append({'path':p.relative_to(ROOT).as_posix(),'before_sha256':sha(backup),'after_sha256':sha(p),'historical_copy':backup.relative_to(ROOT).as_posix()})
    put(HERE/'changes.json',{'scope':'English display labels, one English-guide quotation and prose spacing only. Mathematical statements, evidence bindings, review dates/status decisions, numeric results and original frozen manifests are unchanged. Prior reviewed bytes are retained in before/. Earlier manifests bind those historical versions; this revision binds current presentation.','changes':changes})
    print({'updated_files':len(changes),'review_entries':68,'mathematical_audit_fields_changed':False})
if __name__=='__main__':main()
