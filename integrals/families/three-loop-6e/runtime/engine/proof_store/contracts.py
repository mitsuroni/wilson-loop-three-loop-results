"""Versioned, closed request and receipt contracts; reuse the existing term grammar."""
import hashlib
import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'interfaces'))
from validate import InterfaceError, require, load, digest, file_hash, validate_bundle, term_type, parameters
from schema import obj, array, enum, HASH, TEXT, schema

VERSION = 'wilson-proof-store/0.7'
IR = schema('bundle')
GOAL = IR['properties']['statements']['items']
DEF = IR['properties']['definitions']['items']
TERM = {'$ref': '#/$defs/term'}
RESOURCE = obj(dict(receipt=HASH, arguments=array(TERM), projection=array(enum('left', 'right'))))
AFFINE = obj(dict(format=enum('wilson-proof-request/0.1'), kind=enum('affine_candidate'), parameters=GOAL['properties']['parameters'], assumptions=array(TERM), arguments=array(TERM, minItems=3, maxItems=3), candidate=TERM, evidence=array(HASH, uniqueItems=True)))
STEP_RESOURCE = obj(dict(source={'oneOf':[obj(dict(kind=enum('local'),statement=TEXT)),obj(dict(kind=enum('receipt'),receipt=HASH))]},arguments=array(TERM),projection=array(enum('left','right'))))
PROOF_STEP = obj(dict(goal=GOAL,method=enum('algebra','elementary','positivity','integral','exact','symmetry','transitivity','rewrite','integral_congruence','local_expansion','conjunction'),resources=array(STEP_RESOURCE,maxItems=16)))
REQUEST = {'$schema': IR['$schema'], '$defs': IR['$defs'], 'oneOf': [
    obj(dict(format=enum('wilson-proof-request/0.1'),kind=enum('color_coordinates'),loop_order=dict(type='integer',minimum=1,maximum=8),expression=dict(type='string',minLength=1,maxLength=2048))),
    AFFINE,
    obj(dict(format=enum('wilson-proof-request/0.1'),kind=enum('prove'),definitions=array(DEF,maxItems=64),steps=array(PROOF_STEP,minItems=1,maxItems=32),evidence=array(HASH,uniqueItems=True))),
    obj(dict(format=enum('wilson-proof-request/0.1'),kind=enum('change_variables'),parameters=GOAL['properties']['parameters'],assumptions=array(TERM),integral=TERM,transformed=TERM,mapping=obj(dict(variable=obj(dict(name=TEXT,type=enum('Real'))),body=TERM)),evidence=array(HASH,uniqueItems=True))),
    obj(dict(format=enum('wilson-proof-request/0.1'),kind=enum('benchmark_claims'),benchmark=enum('2e-02'),contract=enum('normalization-v2'),resources=array(HASH,minItems=3,maxItems=3,uniqueItems=True))),
    obj(dict(format=enum('wilson-proof-request/0.1'), kind=enum('integral_equality'), parameters=GOAL['properties']['parameters'], assumptions=array(TERM), integral=TERM, candidate=TERM, evidence=array(HASH,uniqueItems=True))),
    obj(dict(format=enum('wilson-proof-request/0.1'), kind=enum('monomial_candidate'), exponent={'type':'integer','minimum':0,'maximum':1024}, candidate=TERM)),
    obj(dict(format=enum('wilson-proof-request/0.1'), kind=enum('replay'), bundle=IR)),
    obj(dict(format=enum('wilson-proof-request/0.1'), kind=enum('derive'), definitions=array(DEF),
             goal=GOAL, rule=enum('complement_of_sum'), arguments=array(TERM, minItems=2, maxItems=2),
             resources=array(RESOURCE, minItems=1, maxItems=1))),
    obj(dict(format=enum('wilson-proof-request/0.1'), kind=enum('claims'), resources=array(HASH, minItems=1, uniqueItems=True))),
    obj(dict(format=enum('wilson-proof-request/0.1'), kind=enum('analysis'), resources=array(HASH, minItems=5, maxItems=5, uniqueItems=True))),
]}
RECEIPT = obj(dict(format=enum('wilson-proof-receipt/0.1'), goal=HASH, request=HASH, profile=HASH,
                   source=HASH, report=HASH, target=TEXT, run=TEXT,
                   proof_dependencies=array(HASH, uniqueItems=True), planning_dependencies=array(HASH, uniqueItems=True),
                   outcome=enum('checked')))


def canonical(data):
    return json.dumps(data, sort_keys=True, ensure_ascii=True, separators=(',', ':'), allow_nan=False).encode('utf-8')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def validate(data, spec):
    errors = list(Draft202012Validator(spec).iter_errors(data))
    require(not errors, 'contract: ' + (errors[0].message if errors else ''))


def request(data):
    validate(data, REQUEST)
    if data['kind'] == 'replay':
        validate_bundle(data['bundle'])
    if data['kind'] == 'derive':
        defs = {}
        for d in data['definitions']:
            require(d['id'] not in defs, 'duplicate definition')
            require(term_type(d['body'], parameters(d['parameters']), defs) == 'Real', 'definition type')
            defs[d['id']] = len(d['parameters'])
        goal = data['goal']
        env = parameters(goal['parameters'])
        require(not goal['depends_on'], 'use explicit receipt resources')
        for p in goal['assumptions'] + [goal['conclusion']]:
            require(term_type(p, env, defs) == 'Prop', 'goal type')
        for arg in data['arguments'] + data['resources'][0]['arguments']:
            require(term_type(arg, env, defs) == 'Real', 'resource argument type')
    if data['kind'] == 'affine_candidate':
        from .affine import check
        check(data)
    if data['kind'] == 'monomial_candidate':
        from .monomial import check
        check(data)
    if data['kind']=='integral_equality':
        from .integral import check
        check(data)
    if data['kind']=='change_variables':
        from .substitution import check
        check(data)
    if data['kind']=='prove':
        from .proof_steps import check
        check(data)
    return data
