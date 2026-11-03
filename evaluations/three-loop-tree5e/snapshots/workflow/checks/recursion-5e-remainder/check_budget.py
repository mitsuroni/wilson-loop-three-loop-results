"""Check the Laurent-budget hypothesis against immutable targets and fresh sums."""
from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
from compare_t3 import read,sha
from compare_t4_opp import old_sum
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha

def main():
    ptr=read(HERE/'current-t4-budget.json');run=ROOT/ptr['run']
    assert read(run/'process.json')['success']
    fresh=ROOT/read(HERE/'current-t4-final.json')['run']
    assert read(fresh/'process.json')['success']
    bindings={};rows=[]
    def value(name):
        p=run/(name+'.txt');bindings[p.relative_to(ROOT).as_posix()]=sha(p)
        return ha.parse(p)
    def check(name,k,a,b):
        delta=ha.add(a,ha.scale(b,-1))
        rows.append(dict(check=name,power=k,residual=ha.serialize(delta),passed=not delta))
    ids=[13,14,15,33,34,35,40,43]
    for k in range(-6,1):
        old=old_sum(ids,k,bindings);correct=value(f'CorrectEight{k}');bug=value(f'OldBudgetEight{k}')
        p=fresh/f'AdjacentEffective{k}.txt';bindings[p.relative_to(ROOT).as_posix()]=sha(p)
        new=ha.parse(p)
        check('correct_budget_equals_fresh',k,correct,new)
        check('truncation_reproduces_old',k,bug,old)
        check('residual_explained',k,ha.add(correct,ha.scale(bug,-1)),ha.add(new,ha.scale(old,-1)))
        for i in [13,14,15,40]:
            check(f'old_budget_{i}_equals_old',k,value(f'OldBudget{i}q{k}'),old_sum([i],k,bindings))
    out=dict(passed=all(r['passed'] for r in rows),formal='not_run',rows=rows,
             source_sha256=bindings,run=ptr['run'],fresh_run=fresh.relative_to(ROOT).as_posix())
    (HERE/'budget-diagnosis.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('Exact budget diagnosis:',out['passed'])
    for r in rows:
        if not r['passed']:print(r['check'],r['power'],'nonzero coordinates',len(r['residual']))
    return out['passed']

if __name__=='__main__':sys.exit(0 if main() else 1)
