"""Separate maximal-symbol obstructions from period-decorated boundary data."""
from pathlib import Path
import argparse
import hashlib
import json
import compare_weight9 as c
import hlog9 as h

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run');args=ap.parse_args();run=ROOT/args.run
    paths=[run/'gauss-weight9-columns.json',run/'coefficients.json',run/'comparison-gauss-weight9.json']
    assert read(paths[2])['passed']
    cols=[h.lift(row['terms']) for row in read(paths[0])]
    data=read(paths[1])
    targets=[h.lift(data['physical'][name][3]) for name in ('6e-20','6e-21')]
    source=[h.lift(data['bare'][name][3]) for name in ('W','Lsum')]
    rows=[]
    for maximum in (0,2,3,4,5,6,7,8,9):
        def project(v):return {key:co for key,co in v.items() if sum(a*b for a,b in zip(key[1],h.ZW))<=maximum}
        columns=[project(v) for v in cols];rank=len(c.span.build(columns))
        augmented=len(c.span.build(columns+[project(v) for v in targets]))
        srcaug=len(c.span.build(columns+[project(v) for v in source]))
        assert augmented==srcaug
        row={'max_period_weight':maximum,'rank':rank,'augmented_rank':augmented,
             'pair_quotient_dimension':augmented-rank}
        rows.append(row);print(row,flush=True)
    report={'passed':True,'formal':'not_run','weight':9,'rows':rows,
            'meaning':'Maximum zero retains only length-nine Hlog words, i.e. the ordinary maximal-weight symbol coordinates; later rows also retain period-decorated shorter words.',
            'limits':'Formal registered period coordinates. Projection ranks do not assert algebraic independence of period numbers.',
            'inputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in paths}}
    dest=run/'period-filtration.json';assert not dest.exists()
    dest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
