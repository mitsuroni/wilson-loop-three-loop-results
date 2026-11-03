"""Freeze and validate original ordered inputs, never a final-result selector."""
import copy
from collections import OrderedDict
import hashlib
import json
from pathlib import Path
from threading import RLock
from interfaces.validate import load, digest
from preprocessing.ordered import integral_ast,arithmetic
from preprocessing.prepare import verify_prepared
import preprocessing.prepare as _preparation

ROOT = Path(__file__).resolve().parents[1]

# Only fully verified, byte-bound snapshots enter this process-local LRU.
# Keeping one version per path prevents retaining a history of large packages.
_PREPARED_SNAPSHOT_LIMIT = 4
_PREPARED_SNAPSHOT_BYTES = 32 * 1024 * 1024
_PREPARED_SNAPSHOTS = OrderedDict()
_PREPARED_SNAPSHOT_LOCK = RLock()


class InputBudgetError(RuntimeError):pass


def _parse_bytes(raw):
    def pairs(items):
        result={}
        for key,value in items:
            if key in result:raise ValueError('duplicate input JSON key')
            result[key]=value
        return result
    return json.loads(raw.decode('utf-8'),object_pairs_hook=pairs)


def validate_j_budget(j):
    variables=j['variables'];factors=j['propagator_polynomials']
    if len(variables)>64 or len(factors)>32:raise InputBudgetError('endpoint inventory exceeds the bounded recognizer')
    names=set(variables)|{'z1','eps'}
    texts=[j['numerator'],*factors.values(),*(x for pair in j['bounds'].values() for x in pair)]
    if any(not isinstance(s,str) or len(s)>4096 for s in texts):raise InputBudgetError('input expression text budget exceeded')
    def degree(node,depth=0):
        if depth>64:raise InputBudgetError('input AST depth budget exceeded')
        if node['kind']=='rational':
            if len(node['numerator'])>128:raise InputBudgetError('integer size budget exceeded')
            return 0,max(1,abs(int(node['numerator'])).bit_length(),int(node['denominator']).bit_length())
        if node['kind']=='var':return 1,1
        args=node['args'];op=node['op'];bounds=[degree(a,depth+1) for a in args]
        ds=[p[0] for p in bounds];bits=[p[1] for p in bounds]
        if op in ('add','sub'):answer=max(ds),max(bits)+1
        elif op=='neg':answer=bounds[0]
        elif op=='mul':answer=sum(ds),sum(bits)
        elif op=='div' and ds[1]==0:answer=ds[0],sum(bits)
        elif op=='rpow' and args[1]['kind']=='rational' and args[1]['denominator']=='1':
            power=int(args[1]['numerator'])
            if not 0<=power<=64:raise InputBudgetError('exponent budget exceeded')
            answer=ds[0]*power,max(1,bits[0]*power)
        else:raise InputBudgetError('expression outside the bounded polynomial recognizer')
        if answer[1]>4096:raise InputBudgetError('symbolic coefficient bit-size budget exceeded')
        return answer
    def free(node):
        if node['kind']=='var':return {node['name']}
        if node['kind']=='rational':return set()
        return set().union(*(free(a) for a in node['args']))
    for text in factors.values():
        parsed=arithmetic(text,names);used=free(parsed)
        if 'eps' in used or len(used-{'z1'})>2:raise InputBudgetError('pure 2ne propagator endpoint budget exceeded')
        if degree(parsed)[0]>3:raise InputBudgetError('propagator polynomial degree budget exceeded')
    numerator=arithmetic(j['numerator'],names)
    if free(numerator)-{'z1'}:raise InputBudgetError('pure 2ne numerators may depend only on z1')
    if degree(numerator)[0]>64:raise InputBudgetError('numerator degree budget exceeded')
    for pair in j['bounds'].values():
        for text in pair:
            if degree(arithmetic(text,names))[0]>1:raise InputBudgetError('bound polynomial degree budget exceeded')


def _preprocessing_code_hashes():
    """Read the same three actual source files bound by verify_prepared."""
    preparation_path = Path(_preparation.__file__).resolve()
    return {path.relative_to(_preparation.ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (preparation_path, preparation_path.with_name('color.py'),
                         preparation_path.with_name('ordered.py'))}


def _capture_prepared_bytes(path):
    """No stat/mtime shortcut: each call captures all manifest-bound bytes."""
    manifest_raw = path.read_bytes()
    manifest = _parse_bytes(manifest_raw)
    packages = []
    for row in manifest['packages']:
        package_path = (path.parent / row['path']).resolve()
        if package_path.parent != path.parent:
            raise ValueError('package path escapes manifest directory')
        raw = package_path.read_bytes()
        actual_hash = hashlib.sha256(raw).hexdigest()
        if actual_hash != row['sha256']:
            raise ValueError('package changed during input capture')
        packages.append((package_path, row, raw, actual_hash))
    code_hashes = _preprocessing_code_hashes()
    manifest_hash = hashlib.sha256(manifest_raw).hexdigest()
    fingerprint = (manifest_hash, tuple((str(p), h) for p, _, _, h in packages),
                   tuple(sorted(code_hashes.items())))
    return {'manifest_raw': manifest_raw, 'manifest': manifest, 'manifest_hash': manifest_hash,
            'packages': packages, 'code_hashes': code_hashes, 'fingerprint': fingerprint,
            'source_bytes': len(manifest_raw) + sum(len(raw) for _, _, raw, _ in packages)}


def _confirm_prepared_capture(path, captured):
    # Repeat the actual content reads after verification/cache lookup. In
    # particular, the bytes parsed below cannot differ from those verified.
    current = _capture_prepared_bytes(path)
    if current['manifest_raw'] != captured['manifest_raw']:
        raise ValueError('manifest changed during input capture')
    if [(p, raw) for p, _, raw, _ in current['packages']] != [
            (p, raw) for p, _, raw, _ in captured['packages']]:
        raise ValueError('package changed during input capture')
    if current['code_hashes'] != captured['code_hashes']:
        raise ValueError('preprocessing code changed during input capture')
    if captured['manifest'].get('code_sha256') != captured['code_hashes']:
        raise ValueError('preprocessing code checksum mismatch')
    if path.read_bytes() != captured['manifest_raw']:
        raise ValueError('manifest changed during input capture')


def _prepared_snapshot(path, requested_hash=None):
    try:
        captured = _capture_prepared_bytes(path)
        if requested_hash is not None and requested_hash != captured['manifest_hash']:
            raise ValueError('manifest digest mismatch')
        with _PREPARED_SNAPSHOT_LOCK:
            cached = _PREPARED_SNAPSHOTS.get(path)
            if cached is not None and (cached['fingerprint'] != captured['fingerprint']
                                       or cached['manifest_raw'] != captured['manifest_raw']):
                del _PREPARED_SNAPSHOTS[path]
                cached = None
        if cached is None:
            if verify_prepared(path.parent) != captured['manifest']:
                raise ValueError('manifest changed during verification')
        _confirm_prepared_capture(path, captured)
        if cached is not None:
            with _PREPARED_SNAPSHOT_LOCK:
                if path in _PREPARED_SNAPSHOTS:
                    _PREPARED_SNAPSHOTS.move_to_end(path)
            return cached
        by_id = {}
        for _, row, raw, _ in captured['packages']:
            package = _parse_bytes(raw)  # rejects duplicate keys at every depth
            for orbit in package['orbits']:
                by_id.setdefault(orbit['orbit_id'], []).append((package, orbit, row))
        snapshot = {'manifest_raw': captured['manifest_raw'], 'manifest_hash': captured['manifest_hash'],
                    'fingerprint': captured['fingerprint'], 'source_bytes': captured['source_bytes'],
                    'by_id': by_id}
        if snapshot['source_bytes'] <= _PREPARED_SNAPSHOT_BYTES:
            with _PREPARED_SNAPSHOT_LOCK:
                _PREPARED_SNAPSHOTS[path] = snapshot
                _PREPARED_SNAPSHOTS.move_to_end(path)
                while (len(_PREPARED_SNAPSHOTS) > _PREPARED_SNAPSHOT_LIMIT
                       or sum(value['source_bytes'] for value in _PREPARED_SNAPSHOTS.values()) > _PREPARED_SNAPSHOT_BYTES):
                    _PREPARED_SNAPSHOTS.popitem(last=False)
        return snapshot
    except Exception:
        # A detected mutation/failure invalidates even an older matching entry.
        with _PREPARED_SNAPSHOT_LOCK:
            _PREPARED_SNAPSHOTS.pop(path, None)
        raise


def load_input(source, *, root=ROOT):
    if source['kind']=='prepared':
        if set(source)-{'kind','manifest','integral_id','sha256'}: raise ValueError('unknown prepared source field')
        path = Path(source.get('manifest','benchmarks/2ne/prepared/manifest.json'))
        path = path if path.is_absolute() else Path(root)/path
        path = path.resolve()
        if path.name != 'manifest.json': raise ValueError('prepared manifest.json required')
        snapshot=_prepared_snapshot(path,source.get('sha256'))
        matches=snapshot['by_id'].get(source.get('integral_id'),[])
        if len(matches)!=1: raise ValueError('unknown or duplicate integral ID')
        package,orbit,row=matches[0]
        validate_j_budget(orbit['J'])
        return {'kind':'prepared','orbit':copy.deepcopy(orbit),'J_ast':integral_ast(orbit['J']),
                'loop_order':package['loop_order'],'provenance':{'manifest_sha256':snapshot['manifest_hash'],
                'package_sha256':row['sha256'],'orbit_sha256':digest(orbit),'integral_id':orbit['orbit_id']}}
    if set(source)-{'kind','J','loop_order','C','topology','endpoint_edges'}: raise ValueError('unknown explicit source field')
    n=source.get('loop_order')
    if type(n) is not int or n<=0: raise ValueError('positive loop_order required')
    validate_j_budget(source['J']);ast=integral_ast(source['J'])
    orbit={'J':copy.deepcopy(source['J'])}
    for field in ('C','topology','endpoint_edges'):
        if field in source: orbit[field]=copy.deepcopy(source[field])
    return {'kind':'explicit','orbit':orbit,'J_ast':ast,'loop_order':n,
            'provenance':{'explicit_input_sha256':digest(source),'root_admission':'research_input'}}
