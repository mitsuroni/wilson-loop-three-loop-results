"""Closed JSON schemas for the Step 4.1 interchange format (not a proof checker)."""
import json
from pathlib import Path


def obj(properties, required=None):
    return dict(type='object', properties=properties,
                required=list(properties) if required is None else required,
                additionalProperties=False)


def array(items, **kw):
    return dict(type='array', items=items, **kw)


def enum(*values):
    return {'enum': list(values)}


NAME = dict(type='string', pattern=r'^[A-Za-z][A-Za-z0-9_-]*$')
TEXT = dict(type='string', minLength=1)
HASH = dict(type='string', pattern='^[0-9a-f]{64}$')
TERM = {'$ref': '#/$defs/term'}
SOURCE = obj(dict(path=TEXT, sha256_utf8_lf=HASH))
PARAM = obj(dict(name=NAME, type=enum('Real')))


def schema(kind):
    term = {'oneOf': [
        obj(dict(kind=enum('var'), name=NAME)),
        obj(dict(kind=enum('rational'), numerator=dict(type='string', pattern=r'^-?(0|[1-9][0-9]*)$'),
                 denominator=dict(type='string', pattern=r'^[1-9][0-9]*$'))),
        obj(dict(kind=enum('op'), op=enum('add', 'sub', 'mul', 'div', 'neg', 'exp', 'log',
                                         'gamma', 'zeta', 'euler_gamma', 'pi', 'rpow', 'eq', 'ne', 'lt', 'le', 'and', 'implies'),
                 args=array(TERM))),
        obj(dict(kind=enum('call'), name=NAME, args=array(TERM))),
        obj(dict(kind=enum('ite'), condition=TERM, then=TERM, otherwise=TERM)),
        obj(dict(kind=enum('interval_integral'), variable=PARAM, lower=TERM, upper=TERM, body=TERM)),
        obj(dict(kind=enum('cubic_expansion'), variable=PARAM, body=TERM, coefficients=array(TERM,minItems=4,maxItems=4), pole_order=dict(type='integer',minimum=0,maximum=8))),
        obj(dict(kind=enum('local_expansion'), variable=PARAM, body=TERM, constant=TERM, linear=TERM, pole_order=dict(type='integer',minimum=0,maximum=8))),
        obj(dict(kind=enum('has_deriv_at'), variable=PARAM, body=TERM, value=TERM, at=TERM)),
        obj(dict(kind=enum('c1_at'), variables=array(NAME, minItems=1, uniqueItems=True), body=TERM)),
    ]}
    definition = obj(dict(id=NAME, parameters=array(PARAM), body=TERM))
    statement = obj(dict(id=NAME, parameters=array(PARAM), assumptions=array(TERM),
                         conclusion=TERM, depends_on=array(NAME, uniqueItems=True)))
    mapping = obj(dict(statement=NAME, claim_ids=array(TEXT, uniqueItems=True),
                       equations=array(dict(type='integer', minimum=1), uniqueItems=True),
                       lean_source=SOURCE, declaration=TEXT,
                       parameter_groups=array(array(NAME, minItems=1), minItems=1),
                       correspondence=enum('candidate_requires_lean_adapter'), scope_note=TEXT))
    bundle = obj(dict(format=enum('wilson-ir/0.1'), assurance=enum('structural_only'),
                      benchmark=TEXT, source=SOURCE, plan_source=SOURCE,
                      definitions=array(definition), statements=array(statement, minItems=1),
                      mappings=array(mapping)))
    artifact = obj(dict(path=TEXT, sha256_bytes=HASH))
    run = obj(dict(exit_code=dict(type='integer'), stdout=artifact, stderr=artifact,
                   reported_axioms=array(TEXT, uniqueItems=True)))
    certificate = obj(dict(
        format=enum('wilson-certificate/0.1'), status=enum('candidate', 'evidence_unassessed'),
        bundle_sha256=HASH, statement_id=NAME, statement_sha256=HASH,
        definitions_sha256=HASH, assumptions_sha256=HASH, mapping_sha256=HASH,
        upstream=array(obj(dict(statement_id=NAME, certificate_sha256=HASH))),
        checker=obj(dict(profile=enum('lean4-no-sorry-v1'),
                         toolchain=SOURCE, dependency_lock=SOURCE,
                         allowed_axioms={'const': ['propext', 'Classical.choice', 'Quot.sound']})),
        proof={'anyOf': [artifact, {'type': 'null'}]},
        run={'anyOf': [run, {'type': 'null'}]}))
    result = bundle if kind == 'bundle' else certificate
    result.update({'$schema': 'https://json-schema.org/draft/2020-12/schema',
                   '$defs': {'term': term}})
    if kind == 'certificate':
        result['allOf'] = [dict(
            **{'if': {'properties': {'status': {'const': 'candidate'}}}},
            then={'properties': {'run': {'type': 'null'}}},
            **{'else': {'properties': {'proof': artifact, 'run': run}}})]
    return result


if __name__ == '__main__':
    for kind in ('bundle', 'certificate'):
        Path(__file__).with_name(kind + '.schema.json').write_text(
            json.dumps(schema(kind), ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
