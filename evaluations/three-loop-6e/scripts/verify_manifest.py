"""Verify the declared portable packet members, without importing the engine."""
import argparse
import hashlib
import json
from pathlib import Path

def verify(root):
    root=root.resolve();m=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    for rel,row in m['files'].items():
        p=(root/rel).resolve()
        if not p.is_relative_to(root):raise ValueError('Escaping manifest path')
        raw=p.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=row['sha256'] or len(raw)!=row['bytes']:raise ValueError('Changed member '+rel)
    return {'status':'passed','declared_members':len(m['files']),'format':m['format'],
            'scope':'Byte inventory only; not a mathematical proof.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args()
    print(json.dumps(verify(a.root)))
