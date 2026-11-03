"""Verify every declared portable P6 byte against its frozen manifest."""
from pathlib import Path
import hashlib,json
PACK=Path(__file__).resolve().parents[1]
def main():
    manifest=json.loads((PACK/'manifest.json').read_text(encoding='utf-8'));files=manifest['files']
    for rel,r in files.items():
        p=(PACK/rel).resolve()
        assert p.is_relative_to(PACK),rel
        assert p.is_file() and p.stat().st_size==r['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],rel
    print(json.dumps({'status':'passed','declared_members':len(files),'formal':'not_run'}))
if __name__=='__main__':main()
