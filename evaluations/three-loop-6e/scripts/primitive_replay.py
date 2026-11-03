"""Portable original-primitive replay adapted from the prior checked provider.
No installed CAS/library file is read or executed; backend identity is metadata.
"""
from __future__ import annotations
from dataclasses import dataclass
from copy import deepcopy
from pathlib import Path, PureWindowsPath
import hashlib
import json
import ast

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE / 'engine'

@dataclass(frozen=True)
class ArchivedConfig:
    executable: str
    hyperint: str
    periods: str
    backend: dict
    timeout_seconds: float = 120.0
    max_words: int = 50000
    max_output_bytes: int = 8000000
    memory_mb: int = 1024
    def fingerprint(self): return deepcopy(self.backend)

def frozen_config(**limits):
    data=json.loads((PACKAGE/'primitives/archive-profile.json').read_text(encoding='utf-8'))
    b=data['backend_fingerprint']
    return ArchivedConfig(b['executable']['path'],b['hyperint']['path'],b['period_table']['path'],b,**limits)

def historical_maple_path(path):
    value=PureWindowsPath(path).as_posix()
    if any(c in value for c in (chr(34),chr(10),chr(13),chr(0))): raise ValueError('Invalid historical Maple path')
    return value

# Original adapter SHA256: 7c6c1fd55a018947d16bee158773b4eab11cd90b928616f7cf19aae30f0c8c02
UNCHANGED_PRIMITIVE_MODULES = ('backends.py', 'ordered_backend.py', 'algebra.py', 'hlog.py', 'series.py')
RESOURCE_FIELDS = {'timeout_seconds', 'timeout_seconds_per_order', 'total_seconds', 'memory_mb', 'max_words'}
SINGLE_FIELDS = {'operation', 'parameters', 'context', 'definition', 'integration_order',
                 'domain_after_substitution', 'substitution', 'backend'}
RESOURCE_UPDATE_MATH_NODES = {
    'backends.py': ('_script', '_validate_parameters', 'parse_hyperint_coefficient', '_maple_path', '_sha', '_HLOG', '_ZETA'),
    'ordered_backend.py': ('_emit_scalar', 'prepare_ordered_constant', '_script',
                           'prepare_ordered_moments', '_moments_script', '_hash'),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read_bound_file(path, expected):
    raw = Path(path).read_bytes()
    require(sha(raw) == expected, 'archive file digest changed: '+str(path))
    return raw


def match_request(stored, current):
    require(stored['operation'] == current['operation'], 'primitive operation mismatch')
    if stored['operation'] == 'centered-single-opposite-moments':
        require(SINGLE_FIELDS <= stored.keys() and SINGLE_FIELDS <= current.keys(), 'incomplete single profile')
        require({k: stored[k] for k in SINGLE_FIELDS} == {k: current[k] for k in SINGLE_FIELDS},
                'single mathematical request/backend mismatch')
    else:
        require({k: v for k, v in stored.items() if k not in RESOURCE_FIELDS}
                == {k: v for k, v in current.items() if k not in RESOURCE_FIELDS},
                'ordered mathematical request/backend mismatch')


def validate_record(record, script):
    require(record.get('status') == 'computed' and record.get('returncode') == 0
            and record.get('stop_reason') is None, 'primitive execution did not complete')
    require(record.get('script_sha256') == sha(script), 'archived script binding mismatch')


def module_nodes(raw):
    nodes = {}
    for index, node in enumerate(ast.parse(raw.decode('utf-8')).body):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            key = node.name
        elif isinstance(node, ast.Assign):
            key = ','.join(t.id for t in node.targets if isinstance(t, ast.Name)) or f'assignment-{index}'
        else:
            key = f'{type(node).__name__}-{index}'
        nodes[key] = ast.dump(node, include_attributes=False)
    return nodes


def check_resource_update(name, archived, current, approved):
    """Allow only an explicitly reviewed digest, with unchanged math/parser AST."""
    audit = {'module': name, 'archived_sha256': sha(archived), 'current_sha256': sha(current),
             'byte_identical': archived == current}
    if archived == current:
        return audit
    require(name in RESOURCE_UPDATE_MATH_NODES and approved.get(name) == sha(current),
            'primitive implementation changed without exact resource-wrapper approval: '+name)
    old_nodes, new_nodes = module_nodes(archived), module_nodes(current)
    preserved = RESOURCE_UPDATE_MATH_NODES[name]
    for node in preserved:
        require(node in old_nodes and node in new_nodes and old_nodes[node] == new_nodes[node],
                'approved resource update changed mathematical/parser node: '+name+':'+node)
    audit.update(approval='explicit exact module SHA256 supplied to validation CLI',
                 unchanged_mathematical_ast_nodes=list(preserved),
                 changed_top_level_nodes=sorted(k for k in set(old_nodes)|set(new_nodes) if old_nodes.get(k) != new_nodes.get(k)))
    return audit


def single_profile(parameters, backend):
    return {'operation': 'centered-single-opposite-moments', 'parameters': list(parameters),
            'context': 'z > 0; expansion at epsilon=0',
            'definition': '-(1+z)*z^(-eps/2)*integral D^(-1+eps)*t^(a eps)*(1-t)^(b eps)*s^(c eps)*(1-s)^(d eps)',
            'integration_order': ['x', 'y'], 'domain_after_substitution': 'positive quadrant',
            'substitution': 't=x/(1+x), s=y/(1+y)', 'backend': backend}


class PrimitiveReplay:
    def __init__(self, source_run, approved_resources=None):
        from wilson_calc import backends, ordered_backend
        self.source_run = Path(source_run).resolve()
        self.files = {}; self.units = {}; self.events = []; self.scalar_runs = 0; self.module_audits = []
        self.config = frozen_config(max_words=100000, memory_mb=2048)
        self.backend = self.config.fingerprint()
        versions = self.read_json(self.source_run/'source-snapshot/versions.json')
        self.source_versions = versions
        import sympy
        require(versions['sympy'] == sympy.__version__, 'SymPy version differs from primitive profile')
        for name in UNCHANGED_PRIMITIVE_MODULES:
            archived = self.read(self.source_run/'source-snapshot'/name)
            expected = versions['source_sha256'][name]
            require(sha(archived) == expected, 'frozen module differs from source profile')
            self.module_audits.append(check_resource_update(name, archived, (ROOT/'wilson_calc'/name).read_bytes(), approved_resources or {}))
        for name, expected in versions['prepared_sha256'].items():
            require(sha((ROOT/'benchmarks/2ne/prepared'/name).read_bytes()) == expected,
                    'default prepared root changed: '+name)
        for run in sorted((self.source_run/'objects/cas-runs').iterdir()):
            request = self.read_json(run/'request.json')
            result = self.read_json(run/'result.json')
            require(result['status'] == 'computed' and result['provenance'] == request,
                    'primitive result/request binding mismatch')
            require(request['backend'] == self.backend, 'CAS executable/library/period fingerprint mismatch')
            operation = request['operation']
            if operation == 'centered-single-opposite-moments':
                match_request(request, single_profile(request['parameters'], self.backend))
                require(result['complete_orders'] == request['orders'], 'single archive order inventory mismatch')
                for q in result['complete_orders']:
                    records = [r for r in result['records'] if r['order'] == q]
                    require(len(records) == 1, 'duplicate or absent execution record')
                    script = backends._script(request['parameters'], q, self.config).encode('utf-8')
                    self.admit(run/f'order-{q}', run, request, result, records[0], script, [q],
                               {q: {'raw': 'raw.txt', 'reflected': 'reflected.txt'}})
            elif operation in ('ordered-opposite-constant', 'ordered-opposite-gap-moments'):
                if operation == 'ordered-opposite-constant':
                    prepared = ordered_backend.prepare_ordered_constant(request['input_J'], loop_order=request['loop_order'])
                    script = ordered_backend._script(prepared, self.config).encode('utf-8')
                    orders = [0]; names = {0: {'raw': 'raw.txt', 'reflected': 'reflected.txt'}}
                else:
                    prepared = ordered_backend.prepare_ordered_moments(request['input_J'],
                                loop_order=request['loop_order'], max_order=request['max_order'])
                    script = ordered_backend._moments_script(prepared, self.config).encode('utf-8')
                    orders = list(range(request['max_order']+1))
                    names = {q: {c: f'q{q}.{c}.txt' for c in ('raw', 'reflected')} for q in orders}
                match_request(request, {**prepared, 'backend': self.backend})
                require(result['complete_orders'] == orders and len(result['records']) == 1, 'ordered archive inventory mismatch')
                self.admit(run, run, request, result, result['records'][0], script, orders, names)
            else:
                raise ValueError('unrecognized archived primitive operation')

    def read(self, path):
        path = Path(path).resolve()
        require(path.is_relative_to(self.source_run), 'primitive path outside frozen run')
        raw = path.read_bytes()
        self.files[str(path)] = {'sha256': sha(raw), 'bytes': len(raw)}
        return raw

    def read_json(self, path):
        return json.loads(self.read(path).decode('utf-8'))

    def admit(self, folder, run, request, result, record, expected_script, orders, names):
        from wilson_calc.backends import parse_hyperint_coefficient
        script = self.read(folder/'compute.mpl')
        validate_record(record, script)
        require(script == expected_script, 'current generated script differs from archived primitive')
        log = self.read(folder/'maple.log').decode('utf-8')
        marker = {'centered-single-opposite-moments': f'WILSON_MOMENT_COMPLETE order={orders[0]}',
                  'ordered-opposite-constant': 'WILSON_ORDERED_CONSTANT_COMPLETE',
                  'ordered-opposite-gap-moments': 'WILSON_ORDERED_MOMENTS_COMPLETE'}[request['operation']]
        require(marker in log and 'FAILED' not in log and 'Error,' not in log, 'archived completion transcript failed')
        raw_files = {}; parsed = {}
        for q in orders:
            parsed[str(q)] = {}
            for channel, name in names[q].items():
                path = folder/name; raw = self.read(path)
                value = parse_hyperint_coefficient(raw.decode('utf-8'), max_words=request['max_words'])
                require(value == result['coefficients'][str(q)][channel], 'raw primitive differs from archived parsed coefficient')
                raw_files[(q, channel)] = str(path)
                parsed[str(q)][channel] = value
        key = sha(script)
        unit = {'request': request, 'orders': orders, 'folder': str(folder), 'batch_directory': str(run),
                'script_sha256': key, 'raw_files': raw_files, 'original_record': {k: v for k, v in record.items() if k != 'coefficient'}}
        if key in self.units:
            require(self.units[key]['raw_files'] == raw_files, 'duplicate script has multiple archived outputs')
        self.units[key] = unit

    def consume(self, request, script, orders):
        from wilson_calc.backends import parse_hyperint_coefficient
        key = sha(script.encode('utf-8'))
        require(key in self.units, 'no exactly matching archived primitive; CAS fallback is forbidden')
        unit = self.units[key]
        match_request(unit['request'], request)
        require(list(orders) == unit['orders'], 'requested coefficient order differs from archived script')
        coefficients = {}
        for q in orders:
            coefficients[str(q)] = {}
            for channel in ('raw', 'reflected'):
                path = unit['raw_files'][(q, channel)]
                raw = read_bound_file(path, self.files[path]['sha256'])
                coefficients[str(q)][channel] = parse_hyperint_coefficient(raw.decode('utf-8'), max_words=100000)
        event = {'execution_mode': 'cached-primitive-replay', 'operation': request['operation'],
                 'orders': list(orders), 'script_sha256': key, 'source_directory': unit['folder'],
                 'source_request_sha256': self.files[str(Path(unit['batch_directory'])/'request.json')]['sha256'],
                 'raw_files': [{'order': q, 'channel': c, 'path': path, **self.files[path]}
                               for (q, c), path in unit['raw_files'].items()]}
        self.events.append(event)
        return {'status': 'computed', 'coefficients': coefficients, 'complete_orders': list(orders),
                'provenance': {**request, 'primitive_replay': event}, 'run_directory': unit['folder'],
                'records': [{'execution_mode': 'cached-primitive-replay', 'source': unit['original_record']}],
                'formal_state': 'not_run'}

    def single(self, parameters, through, *, config=None, orders=None, **unused):
        from wilson_calc.backends import _script
        config = config or self.config
        request = single_profile(parameters, config.fingerprint())
        selected = list(range(through+1)) if orders is None else list(orders)
        combined = {}; records = []
        for q in selected:
            result = self.consume(request, _script(parameters, q, config), [q])
            combined.update(result['coefficients']); records.extend(result['records'])
        return {'status': 'computed', 'coefficients': combined, 'complete_orders': selected,
                'run_directory': str(self.source_run/'objects/cas-runs'), 'records': records,
                'execution_mode': 'cached-primitive-replay', 'formal_state': 'not_run'}

    def constant(self, j, *, config=None, run_dir=None, cas_budget=None, **kwargs):
        from wilson_calc.ordered_backend import prepare_ordered_constant, _script
        config = config or self.config
        prepared = prepare_ordered_constant(j, **kwargs)
        return self.consume({**prepared, 'backend': config.fingerprint()}, _script(prepared, config), [0])

    def ordered(self, j, *, max_order=2, config=None, run_dir=None, cas_budget=None, **kwargs):
        from wilson_calc.ordered_backend import prepare_ordered_moments, _moments_script
        config = config or self.config
        prepared = prepare_ordered_moments(j, max_order=max_order, **kwargs)
        return self.consume({**prepared, 'backend': config.fingerprint()}, _moments_script(prepared, config), list(range(max_order+1)))

    def report(self):
        return {'format': 'wilson-cached-primitive-replay/0.1', 'source_run': str(self.source_run),
                'source_versions': self.source_versions,
                'verified_unchanged_modules': [r['module'] for r in self.module_audits if r['byte_identical']],
                'primitive_module_audit': self.module_audits,
                'backend_fingerprint': self.backend, 'source_files': self.files,
                'admitted_script_count': len(self.units), 'replay_events': self.events,
                'new_cas_processes': 0, 'formal_state': 'not_run'}

