"""Verify every declared P5 member, without importing mathematical evaluators."""
from pathlib import Path
import hashlib,json
def main():
    root=Path(__file__).resolve().parents[1];m=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    for rel,row in m['files'].items():
        p=(root/rel).resolve()
        if not p.is_relative_to(root):raise ValueError('Escaping declared path')
        raw=p.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=row['sha256'] or len(raw)!=row['bytes']:raise ValueError('Changed declared member '+rel)
    print(json.dumps({'status':'passed','members':len(m['files']),'scope':'Byte inventory, not mathematical proof.'}))
if __name__=='__main__':main()
