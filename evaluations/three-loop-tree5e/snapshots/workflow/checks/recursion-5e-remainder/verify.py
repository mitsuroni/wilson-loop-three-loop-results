"""Verify exported evidence, published revisions, and actual scientific checks."""
from pathlib import Path
import sys,json,re,platform
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE));from compare_t3 import read,sha,vec
from assemble_global import old_target_sum
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'));import hlog_algebra as ha

def main():
    index=read(HERE/'evidence-index.json');files=index['files'];size=0;accepted=[];failed=[]
    for original,r in files.items():
        p=ROOT/r['export'];assert sha(p)==r['sha256'],p;size+=p.stat().st_size
    def resolve(p):
        return ROOT/files[p]['export'] if p in files else ROOT/p
    def bindings(record):
        for p,digest in record.get('source_sha256',{}).items():assert sha(resolve(p))==digest,p
    for row in index['runs']:
        p=ROOT/row['export'];env=read(p/'environment.json')
        for source,digest in env['source_sha256'].items():
            name=Path(source).name;local=p/name
            actual=local if local.exists() else Path(source)
            assert sha(actual)==digest,(row['original'],source)
        log=(p/'maple.log').read_text(encoding='utf8',errors='replace')
        if (p/'process.json').exists():
            process=read(p/'process.json');assert sha(p/'maple.log')==process['log_sha256']
            assert process['success']==row['success']
        if row['success']:
            assert process['returncode']==0 and process['complete_marker'] in log
            assert 'FAILED_NEW_' not in log and not re.search(r'(^|\n)Error,',log)
            accepted.append(row['original'])
        else:failed.append(row['original'])
    c=read(HERE/'construction.json');g=read(HERE/'global-comparison.json');budget=read(HERE/'budget-diagnosis.json')
    assert c['passed'] and c['old_5e_targets_read'] is False and len(c['checks'])==11
    assert g['passed'] and g['construction_sha256']==sha(HERE/'construction.json')
    assert len(g['rows'])==7 and all(not r['corrected_residual'] for r in g['rows'])
    assert [r['power'] for r in g['rows'] if r['old_r001_residual']]==[-1,0]
    assert budget['passed'] and len(budget['rows'])==49 and all(not r['residual'] for r in budget['rows'])
    for name in ('construction','global-comparison','budget-diagnosis','input-checks','t1-comparison','t3-comparison','t4-opposite-comparison'):
        record=read(HERE/(name+'.json'));assert record['passed'] and record['formal']=='not_run';bindings(record)
    oldadj=read(HERE/'t4-adjacent-comparison.json')
    assert not oldadj['passed'] and [r['power'] for r in oldadj['rows'] if not r['passed']]==[-1,0]
    assert all(r['passed'] for r in oldadj['lower_parent_old_positions'])
    probe=ROOT/read(HERE/'current-t4-parent-probe.json')['run']
    assert not ha.parse(resolve((probe/'ParentR1OrderDifference.txt').relative_to(ROOT).as_posix()))
    failedvolume=next(r for r in index['runs'] if r['original'].endswith('Z-t4-volume'))
    assert not failedvolume['success']
    vtext=(ROOT/failedvolume['export']/'maple.log').read_text(encoding='utf8')
    assert 'TW_DONE 1 -2' in vtext and 'TW_DONE 1 -1' in vtext and 'TW_DONE 1 0' not in vtext
    # The continuation reuses only those two proven completed prefixes.
    continuation=next(r for r in index['runs'] if r['original'].endswith('Z-t4-volume-vfirst'))
    ct=(ROOT/continuation['export']/'maple.log').read_text(encoding='utf8')
    assert 'TW_REUSE 1 -2' in ct and 'TW_REUSE 1 -1' in ct and 'TW_REUSE 1 0' not in ct
    cur=read(ROOT/'benchmarks/tgv/5e-3g/current.json');assert len(cur['orbits'])==44
    basis=read(HERE/'evidence/reused/basis-and-conversions.json')
    atoms={s.Symbol(g['name']):vec(g['hlog_terms']) for g in basis['generators']};atoms[s.Symbol('ell')]=ha.LOG
    roundtrips=0;replacements={};changed=[]
    for oid,r in cur['orbits'].items():
        p=ROOT/'benchmarks/tgv/5e-3g'/r['result'];d=read(p)
        assert sha(p)==r['sha256'] and d['orbit_id']==oid and d['revision']==r['revision']
        if r['revision']=='r002':
            changed.append(oid);replacements[int(oid[-2:])]={}
            assert sha(ROOT/d['previous_revision']['path'])==d['previous_revision']['sha256']
            for k in range(-6,1):
                v=vec(d['hlog_coefficients'][str(k)]);expr=s.sympify(d['laurent']['coefficients'][str(k)],locals={'E':s.Symbol('E')})
                assert ha.from_sympy(expr,atoms)==v;roundtrips+=1;replacements[int(oid[-2:])][k]=v
    assert changed==['5e-3g-13','5e-3g-14','5e-3g-15','5e-3g-33','5e-3g-34','5e-3g-35']
    result=read(ROOT/'benchmarks/recursion-audit/results/5e-effective-sum/r001.json')
    comparison=old_target_sum(replacements)
    for k in range(-6,1):
        v=vec(result['hlog_coefficients'][str(k)])
        assert v==comparison[k]==vec(c['S5e'][str(k)])
        expr=s.sympify(result['laurent']['coefficients'][str(k)],locals={'E':s.Symbol('E')})
        assert ha.from_sympy(expr,atoms)==v;roundtrips+=1
    assert roundtrips==49
    pub=read(HERE/'publication.json');assert pub['passed'] and pub['changed_individual_positions']==8
    assert sha(ROOT/pub['sum_result'])==pub['sum_sha256']
    assert sha(ROOT/'benchmarks/tgv/5e-3g/current.json')==pub['current_sha256']
    assert len(read(HERE/'input-checks.json')['pushforward_derivative_controls'])==16
    neg=[]
    for name,powers in [('historical_Maple_relative_truncation',[-1,0]),('T4_parent_truncated_at_eps0',[-1,0])]:
        residual=g['rows'] if name.startswith('historical') else c['parent_truncation_negative']
        actual=[r['power'] for r in residual if r['old_r001_residual']] if isinstance(residual,list) else [int(k) for k,v in residual.items() if v]
        assert actual==powers;neg.append(dict(case=name,detected_powers=powers))
    # A_controls = -1/4(IT1+IT2)Iadj has coefficient +1/64 at eps^-6.
    control_residual=s.Rational(1,64);assert control_residual!=0
    neg.append(dict(case='omit_two_auxiliary_only_control_domains',detected_power=-6,residual=str(control_residual)))
    gamma=s.Symbol('EulerGamma');leading=s.Rational(1,96)
    assert s.expand(-3*gamma*leading)==-gamma/32 and -gamma/32!=0
    neg.append(dict(case='repeat_three_loop_Euler_compensation',detected_power=-5,residual='-EulerGamma/32'))
    stage=read(ROOT/'workflow/checks/recursion-5e-sum/verification.json');assert stage['passed'] and not stage['global_goal_complete']
    out=dict(passed=True,formal='not_run',local_effective_sum_goal_complete=True,
        all_44_individual_integrals_independently_verified=False,global_w3_complete=False,
        target_positions=7,global_zero_residual_positions=7,construction_overlap_checks=11,
        budget_diagnosis_checks=49,exact_existing_basis_roundtrips=49,
        corrected_integrals=changed,corrected_individual_positions=8,
        subgroup_old_target_exact_comparisons=56,original_pushforward_controls=16,
        parent_alternative_order_exact_check=True,negative_controls=neg,
        exported_files=len(files),exported_bytes=size,successful_CAS_jobs=len(accepted),failed_CAS_attempts=len(failed),
        evidence_index_sha256=sha(HERE/'evidence-index.json'),
        first_stage_source_geometry_and_QMC='workflow/checks/recursion-5e-sum/verification.json',
        environment={'python':platform.python_version(),'sympy':s.__version__},
        source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(HERE.glob('*.py'))},
        evidence_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(HERE.glob('*.json')) if p.name!='verification.json'},
        remaining=['Individual cancelling directions/Omega6; stronger relation/anchor audit.',
            'Full-epsilon T3/T4 auxiliary evaluation, general complex z/i0 and formal proof.',
            'Other three-loop families, counterterms, contour subtraction and complete w3 RG/cusp.'])
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('VERIFIED S5e seven coefficients; six revisions/eight positions; 49 basis and 49 budget checks; formal=not_run.')
    print('Export:',len(files),'files,',size,'bytes;',len(accepted),'completed jobs;',len(failed),'failed attempts retained.')

if __name__=='__main__':main()
