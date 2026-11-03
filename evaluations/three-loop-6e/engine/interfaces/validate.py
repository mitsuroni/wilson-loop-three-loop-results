"""Validate syntax, types, bindings and fingerprints. NEVER accept a mathematical claim."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


class InterfaceError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InterfaceError(message)


def load(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key: ' + key)
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs)


def digest(value):
    # v0.1 canonical JSON: ordered arrays, sorted ASCII keys, no floats in accepted IR.
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=True,
                                    separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest()


def file_hash(path, text=True):
    payload = Path(path).read_text(encoding='utf-8').replace('\r\n', '\n').encode('utf-8') if text else Path(path).read_bytes()
    return hashlib.sha256(payload).hexdigest()


def source_path(source, root):
    name = source['path']
    require('\\' not in name and ':' not in name and not name.startswith('/'), 'path must be repository-relative POSIX')
    path = (root / name).resolve()
    require(path.is_relative_to(root.resolve()), 'path escapes repository')
    require(path.is_file(), 'missing file: ' + name)
    return path


def check_source(source, root):
    key = 'sha256_utf8_lf' if 'sha256_utf8_lf' in source else 'sha256_bytes'
    require(file_hash(source_path(source, root), key == 'sha256_utf8_lf') == source[key],
            'stale source: ' + source['path'])


def check_schema(data, kind):
    spec = load(HERE / (kind + '.schema.json'))
    Draft202012Validator.check_schema(spec)
    errors = list(Draft202012Validator(spec).iter_errors(data))
    require(not errors, 'schema: ' + (errors[0].message if errors else ''))


def parameters(items):
    names = [p['name'] for p in items]
    require(len(names) == len(set(names)), 'duplicate parameter')
    return set(names)


def term_type(node, env, definitions):
    kind = node['kind']
    if kind == 'var':
        require(node['name'] in env, 'unbound variable: ' + node['name'])
        return 'Real'
    if kind == 'rational':
        n, d = int(node['numerator']), int(node['denominator'])
        require(str(n) == node['numerator'] and math.gcd(n, d) == 1, 'noncanonical rational')
        return 'Real'
    if kind in ('op', 'call'):
        args = [term_type(a, env, definitions) for a in node['args']]
        if kind == 'call':
            require(node['name'] in definitions, 'unknown or forward definition: ' + node['name'])
            require(args == ['Real'] * definitions[node['name']], 'definition arity/type mismatch')
            return 'Real'
        op = node['op']
        arity = 0 if op in ('pi', 'euler_gamma') else (1 if op in ('neg', 'exp', 'log', 'gamma', 'zeta') else 2)
        expected = 'Prop' if op in ('and', 'implies') else 'Real'
        require(args == [expected] * arity, 'operator arity/type mismatch: ' + op)
        return 'Prop' if op in ('eq', 'ne', 'lt', 'le', 'and', 'implies') else 'Real'
    if kind == 'ite':
        require(term_type(node['condition'], env, definitions) == 'Prop', 'ite condition must be Prop')
        for key in ('then', 'otherwise'):
            require(term_type(node[key], env, definitions) == 'Real', 'ite branch must be Real')
        return 'Real'
    if kind == 'c1_at':
        require(set(node['variables']) <= env, 'unknown C1 coordinate')
        require(term_type(node['body'], env, definitions) == 'Real', 'C1 body must be Real')
        return 'Prop'
    variable = node['variable']['name']
    require(variable not in env, 'binder shadows existing variable')
    require(term_type(node['body'], env | {variable}, definitions) == 'Real', 'binder body must be Real')
    if kind == 'cubic_expansion':
        for c in node['coefficients']:
            require(term_type(c, env, definitions) == 'Real', 'coefficient must be Real')
        return 'Prop'
    keys = ('constant', 'linear') if kind == 'local_expansion' else (('lower', 'upper') if kind == 'interval_integral' else ('value', 'at'))
    for key in keys:
        require(term_type(node[key], env, definitions) == 'Real', 'bound/value must be Real')
    return 'Real' if kind == 'interval_integral' else 'Prop'


def validate_bundle(data, root=ROOT):
    check_schema(data, 'bundle')
    check_source(data['source'], root)
    check_source(data['plan_source'], root)
    plan = load(source_path(data['plan_source'], root))
    require(plan['benchmark_id'] == data['benchmark'], 'benchmark mismatch')
    require(plan['source']['sha256_utf8_lf'] == data['source']['sha256_utf8_lf'], 'plan derivation mismatch')
    claims = {c['id'] for c in plan['claims']}
    definitions, statements = {}, {}
    for definition in data['definitions']:
        name = definition['id']
        require(name not in definitions, 'duplicate definition')
        env = parameters(definition['parameters'])
        require(term_type(definition['body'], env, definitions) == 'Real', 'definition must return Real')
        definitions[name] = len(env)
    for statement in data['statements']:
        name = statement['id']
        require(name not in statements, 'duplicate statement')
        env = parameters(statement['parameters'])
        for prop in statement['assumptions'] + [statement['conclusion']]:
            require(term_type(prop, env, definitions) == 'Prop', 'statement requires Prop')
        require(set(statement['depends_on']) <= statements.keys(), 'dependency must precede statement (no cycles)')
        statements[name] = statement
    mapped = set()
    for mapping in data['mappings']:
        name = mapping['statement']
        require(name in statements and name not in mapped, 'unknown or duplicate mapping')
        mapped.add(name)
        require(set(mapping['claim_ids']) <= claims, 'unknown benchmark claim')
        require(set(mapping['equations']) <= set(plan['source']['equations']), 'unknown source equation')
        flat = [n for group in mapping['parameter_groups'] for n in group]
        require(flat == [p['name'] for p in statements[name]['parameters']], 'Lean parameter layout mismatch')
        check_source(mapping['lean_source'], root)
    require(mapped == statements.keys(), 'each statement needs an explicit candidate mapping')
    return statements


def bindings(bundle, statement_id):
    statement = next(s for s in bundle['statements'] if s['id'] == statement_id)
    mapping = next(m for m in bundle['mappings'] if m['statement'] == statement_id)
    return dict(bundle_sha256=digest(bundle), statement_id=statement_id,
                statement_sha256=digest(statement), definitions_sha256=digest(bundle['definitions']),
                assumptions_sha256=digest(statement['assumptions']), mapping_sha256=digest(mapping))


def validate_certificate(data, bundle, root=ROOT):
    check_schema(data, 'certificate')
    statements = validate_bundle(bundle, root)
    require(data['statement_id'] in statements, 'unknown certificate statement')
    for key, value in bindings(bundle, data['statement_id']).items():
        require(data[key] == value, 'certificate binding mismatch: ' + key)
    upstream = [u['statement_id'] for u in data['upstream']]
    require(upstream == statements[data['statement_id']]['depends_on'], 'upstream dependency mismatch')
    for key in ('toolchain', 'dependency_lock'):
        check_source(data['checker'][key], root)
    if data['proof'] is not None:
        check_source(data['proof'], root)
    if data['run'] is not None:
        for key in ('stdout', 'stderr'):
            check_source(data['run'][key], root)
    # Exit code / axioms / upstream hashes are untrusted evidence at this stage.
    # There is deliberately no transition to accepted, even for exit_code == 0.
    return 'structurally_valid_not_verified'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundle', type=Path)
    parser.add_argument('--certificate', type=Path)
    args = parser.parse_args()
    try:
        bundle = load(args.bundle)
        statements = validate_bundle(bundle)
        if args.certificate:
            validate_certificate(load(args.certificate), bundle)
        print(json.dumps(dict(status='structurally_valid_not_verified', statements=len(statements), accepted_claims=[])))
    except (InterfaceError, OSError, ValueError, KeyError) as exc:
        parser.exit(1, str(exc) + '\n')


if __name__ == '__main__':
    main()
