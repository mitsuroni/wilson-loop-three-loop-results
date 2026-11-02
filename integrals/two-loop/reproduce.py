"""Rebuild two-loop series and exact representations into a fresh directory."""
from pathlib import Path
import argparse,hashlib,importlib.util,json,re,sys,time
from fractions import Fraction

sys.dont_write_bytecode=True
import sympy as s
from sympy.parsing.mathematica import parse_mathematica
from family_engines import reconstruct,canonical,parse,eps,ell,M,N

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def dump(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    out=a.output_dir.resolve()
    public_root=HERE.parents[1]
    if out.is_relative_to(public_root) and not out.is_relative_to(public_root/'replay-output'):
        raise SystemExit('Inside this repository choose replay-output/<new-name>; immutable science/docs/package trees cannot be output directories.')
    if out.exists() and any(out.iterdir()):raise SystemExit('Choose a new or empty output directory; frozen science is not overwritten.')
    out.mkdir(parents=True,exist_ok=True);tick=time.monotonic();checks=[]
    manifest=read(HERE/'source-manifest.json')
    for rel,v in manifest['sources'].items():assert sha(HERE/rel)==v['sha256'],rel
    spec=importlib.util.spec_from_file_location('historical_hlog_algebra',HERE/'sources/workflow/checks/single-opposite/hlog_algebra.py')
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    atoms={s.Symbol(k):h.parse(HERE/v) for k,v in manifest['moment_files'].items()}
    basis_data=read(HERE/'sources/benchmarks/2ne/families/single-opposite/basis-functions.json')
    basis=basis_data.get('functions',basis_data)['tauDq4']
    tau={(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in basis['hlog_terms']}
    atoms[s.Symbol('tau_D4')]=tau;atoms[ell]=h.LOG
    def vector(x):
        x=s.expand(x)
        for sym in list(x.free_symbols):
            m=re.fullmatch(r'N_(\d+)_(\d+)',str(sym))
            if m:x=x.subs(sym,M(int(m[1]),int(m[2])+1)-ell*M(int(m[1]),int(m[2])))
        return h.from_sympy(s.expand(x),atoms)
    def eq(name,v):
        ok=not v;checks.append({'name':name,'passed':ok});assert ok,name
    compact=read(HERE/'sources/workflow/checks/tgv04-high-order/compact.json')
    cs=reconstruct(compact)
    refs={}
    p=read(HERE/'sources/workflow/checks/two-loop-report-20261015/pure-data.json')
    refs.update({oid:{int(q):parse(v['sympy']) for q,v in row.items()} for oid,row in p['coefficients'].items()})
    t=read(HERE/'sources/workflow/checks/two-loop-report-20261015/tgv-data.json')
    refs.update({oid:{int(q):parse(v['expression']) for q,v in row.items()} for oid,row in t['integrals'].items()})
    b=read(HERE/'sources/workflow/checks/two-loop-report-20261015/self-energy-data.json')
    refs.update({oid:{int(q):parse(v['sympy']) for q,v in row.items()} for oid,row in b['scalar_coefficients'].items()})
    results={};rows=[]
    for oid,row in cs.items():
        source=manifest['canonical_sources'][oid];old=read(HERE/source['path'])
        for q in range(-4,3):eq(f'{oid} rebuilt mother coefficient versus previous report eps{q}',h.add(vector(row[q]),h.scale(vector(refs[oid][q]),-1)))
        for q in range(-4,1):
            entry=old['laurent']['coefficients'][str(q)]
            if 'hlog_terms' in entry:
                target={(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in entry['hlog_terms']}
            else:
                expression=parse_mathematica(entry['expression']) if entry['language'].startswith('Wolfram') else parse(entry['expression'])
                expression=expression.subs({s.Symbol('Pi'):s.pi,s.log(s.Symbol('z')):ell})
                expression=expression.xreplace({v:s.zeta(v.args[0]) for v in expression.atoms(s.Function) if str(v.func)=='Zeta'})
                target=vector(expression)
            eq(f'{oid} current individual source eps{q}',h.add(vector(row[q]),h.scale(target,-1)))
        coeffs={}
        for q,v in row.items():
            # Use the old compact representative only after full exact equality;
            # construction above never uses the pure/Pi1 coefficient table.
            published=refs[oid][q]
            coeffs[str(q)]={'sympy':str(published),'latex':latex(published),'hlog_terms':h.serialize(vector(v))}
        fam='pure' if oid.startswith('4e-') else 'TGV' if oid.startswith('3e-') else 'Pi1'
        result={'format':'wilson-two-loop-integral/1','orbit_id':oid,'family':fam,'domain':'fixed real z>0; inherited source i0; no new full complex-domain theorem','canonical_convention':source['canonical_convention'],'assembly_sign':1,'legacy_to_canonical_multiplier':-1 if fam=='Pi1' else 1,'color_weights_over_CF':{'CA':'(10-6*eps)/(3-2*eps)','T':'-8*(1-eps)/(3-2*eps)'} if fam=='Pi1' else {'CA':'1','T':'0'},'coefficients':coeffs,'powers':list(range(-4,3)),'remainder':'O(eps^3)','canonical_source':source,'calculation_this_run':'Evaluated mother-function Taylor/Gamma convolution plus exact whole-HLOG comparison; source integration proof and HLOG moment dictionary reused' if oid!='3e-3g-04' else 'Historical complete HyperInt-integrated coefficients are reused; this run verifies full HLOG-to-MN/tau representations and source bindings. No new Maple integral evaluation.','formal':'not_run'}
        dump(out/oid/'result.json',result);results[oid]=result
        rows.append({k:result[k] for k in ['orbit_id','family','canonical_convention','assembly_sign','legacy_to_canonical_multiplier','color_weights_over_CF','canonical_source']}|{'result':oid+'/result.json'})
    for q in (1,2):
        integrated=h.parse(HERE/f'sources/workflow/checks/tgv04-high-order/runs/positive-epsilon/I{q}.txt')
        eq(f'3e-3g-04 complete historical integrated HLOG eps{q}',h.add(vector(cs['3e-3g-04'][q]),h.scale(integrated,-1)))
    # Fresh finite-order Q-to-I assembly from integrated Q0,Q1,Q2.
    qv={q:h.parse(HERE/manifest['four_e04_Q_files'][str(q)]) for q in range(3)}
    p2=s.zeta(2)-ell**2/2
    for q,target in [(0,h.scale(qv[0],Fraction(-1,8))),(1,h.scale(qv[1],Fraction(-1,8))),(2,h.scale(h.add(qv[2],h.mul(vector(p2),qv[0])),Fraction(-1,8)))]:
        eq(f'4e-04 integrated two-variable Q to physical I eps{q}',h.add(vector(cs['4e-04'][q]),h.scale(target,-1)))
    # Wrong signs/multiplicities must fail in actual complete coordinates.
    assert vector(2*cs['4e-04'][0])!=vector(cs['4e-04'][0])
    assert vector(-cs['2e-Pi1-01'][-3])!=vector(cs['2e-Pi1-01'][-3])
    summary={'status':'passed','checks':len(checks),'negative_controls':2,'integrals':12,'coefficients':84,'source_files':len(manifest['sources']),'seconds':time.monotonic()-tick,'checks_detail':checks,'new_original_12_graph_integrations':False,'new_Maple_or_Mathematica_execution':False,'source_integration_and_HLOG_dictionary_reused':True,'formal':'not_run'}
    dump(out/'catalog.json',{'format':'wilson-two-loop-calculation-catalogue/1','path_base':'integrals/two-loop','expression_language':'SymPy exact strings, eps; raw M_k_j/N_k_j moment symbols and zeta(n); tau_D4 at TGV04 eps2','external_CF_removed':True,'canonical_convention':'All I include one Euler compensation and recorded C4 normalization. Native Pi1 r002 already carries the physical sign: consume +B1, never reapply legacy -1.','integrals':rows,'remainder':'O(eps^3)','validation':summary})
    dump(out/'verification.json',summary)
    print(json.dumps({k:summary[k] for k in ['status','checks','negative_controls','integrals','coefficients','seconds','formal']}))

def latex(x):
    names={a:(r'\tau^+_{D,4}' if str(a)=='tau_D4' else re.sub(r'([MN])_(\d+)_(\d+)',r'\1_{\2,\3}',str(a))) for a in x.free_symbols}
    return s.latex(x,symbol_names=names)

if __name__=='__main__':main()
