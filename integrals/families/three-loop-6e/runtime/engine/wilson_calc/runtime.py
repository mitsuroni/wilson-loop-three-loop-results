"""Content-addressed candidate calculations; no proof database or receipts."""
import hashlib
import json
from pathlib import Path
import tempfile
import os
import subprocess
import sys
from uuid import uuid4
import sympy
from interfaces.validate import digest, load
from proof_store.storage import Store

ROOT = Path(__file__).resolve().parents[1]


def scalar_worker(payload,seconds,memory_mb):
    """Bound expensive symbolic expansion in a separate process."""
    from .backends import run_bounded_process
    folder=ROOT/'.wilson/calculations'/('scalar-'+uuid4().hex);folder.mkdir(parents=True)
    source=folder/'request.json';source.write_text(json.dumps(payload),encoding='utf-8')
    log=folder/'worker.log'
    process=run_bounded_process([sys.executable,'-X','utf8','-m','wilson_calc._worker',str(source)],cwd=ROOT,
                  stdout_path=log,timeout_seconds=seconds,memory_mb=memory_mb)
    if process['stop_reason'] in ('time','memory'):raise TimeoutError('scalar worker exceeded '+process['stop_reason']+' budget')
    stdout=log.read_text(encoding='utf-8')
    if process['returncode']:raise ArithmeticError('scalar worker failed: '+stdout[-1000:])
    result=json.loads(stdout)
    if result.get('format')!='wilson-scalar-worker/0.1':raise ValueError('invalid scalar worker response')
    return result['coefficients'],result['lower_bound']


class CalculationObjects(Store):
    """Use existing immutable object methods without opening its receipt index."""
    def __init__(self,path):
        self.path=Path(path).resolve()
        if not self.path.is_relative_to(ROOT): raise ValueError('calculation store must remain in workspace')
        (self.path/'objects').mkdir(parents=True,exist_ok=True)

    def close(self): pass

    def put(self,data):
        try:return super().put(data)
        except PermissionError:
            if sys.platform!='win32':raise
        # Windows restricted filesystems may allow atomic rename but prohibit
        # hard links. Windows rename never replaces an existing destination.
        key=hashlib.sha256(data).hexdigest();path=self.object_path(key)
        path.parent.mkdir(exist_ok=True)
        if path.exists():
            if self.read(key)!=data:raise ValueError('candidate object collision')
            return key
        fd,name=tempfile.mkstemp(dir=path.parent,prefix='.candidate-')
        try:
            with os.fdopen(fd,'wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
            try:os.rename(name,path)
            except FileExistsError:
                if self.read(key)!=data:raise ValueError('candidate object collision')
        finally:Path(name).unlink(missing_ok=True)
        return key


def engine_fingerprint():
    paths=sorted(Path(__file__).parent.glob('*.py'))
    return digest({'source':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if not p.name.startswith('test_')},
                   'sympy':sympy.__version__,'hlog_kernel':'dt/(t-a)','endpoint_convention':'ordinary-or-explicit-tangential/0.1'})


class Runtime:
    def __init__(self,store=None):
        self.store=store or CalculationObjects(ROOT/'.wilson/calculation-cache')
        self.cache=self.store.path/'lookup';self.cache.mkdir(exist_ok=True)

    def read(self,key):
        return self.store.json(key)

    def lookup(self,semantic,through,quantity):
        index=self.cache/(semantic+'.json')
        if not index.exists():return None
        row=load(index); result=self.read(row['digest'])
        from .contracts import validate_result
        validate_result(result,semantic_key=semantic,through=through,quantity=quantity)
        if result.get('status')=='computed' and result.get('through')==through:return result
        return None

    def save(self,result):
        # Save checkpoint body before adding its own digest; avoid a circular hash.
        body=json.loads(json.dumps(result));body.get('provenance',{}).pop('checkpoint',None)
        key=self.store.put_json(body)
        if result['status']=='computed':
            destination=self.cache/(result['provenance']['semantic_key']+'.json')
            fd,name=tempfile.mkstemp(dir=self.cache,prefix='.lookup-')
            try:
                with os.fdopen(fd,'w',encoding='utf-8') as f:json.dump({'digest':key},f)
                os.replace(name,destination)
            finally:Path(name).unlink(missing_ok=True)
        result['provenance']['checkpoint']=key
        return result
