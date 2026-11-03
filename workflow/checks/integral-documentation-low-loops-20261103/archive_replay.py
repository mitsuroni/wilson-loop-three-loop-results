"""Preserve compact actual outputs and compare primary and copied executions."""
import shutil
from prepare import HERE,ROOT,READING,read,put,sha

def main():
    out=HERE/'evidence/replay';out.mkdir(parents=True,exist_ok=True)
    iso=HERE/'runs/isolation-r001';assert read(iso/'replay.json')['status']=='passed'
    for family,names in [('one-loop',['checks.json','numeric.json','generated-before-holdout.json']),('two-loop',['verification.json','catalog.json']),('two-loop-controls',['integral-checks.json'])]:
        for name in names:
            p=iso/family/name;dst=out/family/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
    for name in ('one-loop','two-loop','two-loop-controls'):
        shutil.copyfile(iso/(name+'.log'),out/(name+'-isolated.log'))
    # The first adapter failure is retained as an implementation finding.
    shutil.copyfile(HERE/'two-loop.log',out/'two-loop-first-import-failure.log')
    for stage in ('one-loop-r001','two-loop-r002','two-loop-controls-r001'):
        for p in (HERE/'runs'/stage).rglob('*'):
            if p.is_file():
                dst=HERE/'evidence/primary'/stage/p.relative_to(HERE/'runs'/stage)
                dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
    members=[]
    for row in read(HERE/'entries.json'):
        family=row['loop_directory'];oid=row['orbit_id']
        new=iso/family/oid/'result.json'
        first=HERE/'runs'/('one-loop-r001' if family=='one-loop' else 'two-loop-r002')/oid/'result.json'
        assert new.read_bytes()==first.read_bytes(),oid
        assert new.read_bytes()==(ROOT/row['path']/'result.json').read_bytes(),oid
        dst=out/family/oid/'result.json';dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(new,dst)
        obj=read(new);members.append({'orbit_id':oid,'sha256':sha(new),'primary_and_copy_same_bytes':True,'adopted_reader_same_bytes':True,'coefficients':len(obj['laurent']['coefficients']) if 'laurent' in obj else len(obj['coefficients'])})
    one=read(iso/'one-loop/checks.json');two=read(iso/'two-loop/verification.json');ctrl=read(iso/'two-loop-controls/integral-checks.json');numeric=read(iso/'one-loop/numeric.json')
    report={'status':'passed','integrals':len(members),'coefficients':sum(r['coefficients'] for r in members),
        'one_loop_checks':one['check_count'],'one_loop_numeric_controls':len(numeric['controls']),'one_loop_negative_controls':one['negative_controls'],
        'two_loop_checks':two['checks'],'two_loop_negative_controls':two['negative_controls'],'two_loop_mathematical_controls':ctrl['check_count'],
        'copied_execution':read(iso/'replay.json'),'members':members,
        'initial_implementation_failure':{'problem':'New runpy adapter lacked sibling-module import path for family_engines.','fix':'Added shared runtime to sys.path in new wrappers only.','retest':'Affected primary two-loop entrypoint and then all three stages from a real copied closure passed.','log':'evidence/replay/two-loop-first-import-failure.log'},
        'scope':'One-loop source regeneration and two-loop Gamma/mother/adopted HLOG assembly with historical controls. No new full TGV04 five-dimensional CAS integration or fourteen independent cold integrations. This report does not grade prose or grant formal evidence.','formal':'not_run'}
    put(HERE/'replay.json',report)
    print({k:report[k] for k in ('status','integrals','coefficients','one_loop_checks','one_loop_numeric_controls','two_loop_checks','two_loop_mathematical_controls')})
if __name__=='__main__':main()
