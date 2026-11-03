"""Point new articles to retained outputs, rather than ignored run directories."""
from prepare import ROOT,READING

def main():
    changed=[]
    for folder in ('one-loop','two-loop'):
        for p in (READING/folder).rglob('*.md'):
            old=p.read_text(encoding='utf-8-sig');new=old
            for stage in ('one-loop-r001','two-loop-r002','two-loop-controls-r001'):
                new=new.replace('20261103/runs/'+stage,'20261103/evidence/primary/'+stage)
            if new!=old:p.write_text(new,encoding='utf-8');changed.append(p.relative_to(ROOT).as_posix())
    print({'article_receipt_links_retained':len(changed),'mathematical_text_changed':False})
if __name__=='__main__':main()
