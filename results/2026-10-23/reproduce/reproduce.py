"""Replay the frozen standalone delivery into a new directory."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
    scripts=Path(__file__).resolve().parent;package=scripts.parent;root=package/'sources';out=args.output_dir.resolve()
    if out.exists() and any(out.iterdir()): raise SystemExit('Choose a new or empty output directory; frozen evidence is not overwritten.')
    out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((package/'source-manifest.json').read_text(encoding='utf-8'))['source_sha256']
    for rel,h in manifest.items():
        assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==h, ('source hash mismatch',rel)
    runs=[]
    def run(script,*options):
        command=[sys.executable,'-X','utf8',str(script),*map(str,options)]
        result=subprocess.run(command,text=True,encoding='utf-8',capture_output=True)
        runs.append({'script':str(script.relative_to(package)), 'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
        print(result.stdout,end='',flush=True)
        if result.returncode: print(result.stderr,file=sys.stderr); raise RuntimeError(script)
    run(root/'workflow/checks/matter-rg-20261017/check.py','--root',root,'--output-dir',out/'rg-replay')
    run(scripts/'build_coefficients.py','--root',root,'--output-dir',out,'--replay-dir',out/'rg-replay')
    run(scripts/'build_diagrams.py','--root',root,'--output-dir',out)
    run(scripts/'verify_uv.py','--root',root,'--output-dir',out)
    run(scripts/'numeric_samples.py','--output-dir',out)
    original=json.loads((package/'coefficients.json').read_text(encoding='utf-8'))
    rebuilt=json.loads((out/'coefficients.json').read_text(encoding='utf-8'))
    assert original['channels']==rebuilt['channels'], 'Rebuilt channel data changed'
    assert original['lower_orders']==rebuilt['lower_orders'], 'Rebuilt lower-order data changed'
    a=json.loads((package/'diagram-ledger.json').read_text(encoding='utf-8'))
    b=json.loads((out/'diagram-ledger.json').read_text(encoding='utf-8'))
    assert a['diagrams']==b['diagrams'], 'Rebuilt diagram data changed'
    (out/'reproduction.json').write_text(json.dumps({'status':'passed','date':'2026-10-23','sources_verified':len(manifest),'channel_and_diagram_semantic_equality':True,'commands':runs},indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','sources':len(manifest),'reproduced_channels':3,'reproduced_diagrams':31}))

if __name__=='__main__': main()
