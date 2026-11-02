"""Physics controls for the common convention, independent of orbit identifiers."""
from copy import deepcopy
import pytest
import sympy as s
from montopar.sign_rules import (expected_contract,validate_contract,coordinate_phase,
                                 three_slot_trace_phase,transport_word)
from montopar.c_factor import compute_c_tgv,compute_c_fgv
from montopar.emitter import _combine_c_tgv,_build_c_fgv_block

def expression(text):
    return s.sympify(text.replace('[','(').replace(']',')').replace('^','**'),
        locals={'I':s.I,'Pi':s.pi,'Gamma':s.gamma})

def source():
    raw={'vertices':[{'line':3,'color':'a1'},{'line':4,'color':'a2'},{'line':4,'color':'a3'}]}
    raw['physical_convention']=expected_contract(raw['vertices'])
    return raw

def test_transport_word_and_actual_SU2_f_trace():
    raw=source();word=transport_word(raw['vertices'])
    assert word==['a3','a2','a1']
    sigma=[s.Matrix([[0,1],[1,0]])/2,s.Matrix([[0,-s.I],[s.I,0]])/2,s.diag(1,-1)/2]
    value=sum(s.LeviCivita(a,b,c)*s.trace(sigma[c]*sigma[b]*sigma[a])/2
              for a in range(3) for b in range(3) for c in range(3))
    assert value==-3*s.I/4
    assert three_slot_trace_phase(raw)==-1
    raw['orbit_id']='arbitrary label; no family sign lookup'
    assert three_slot_trace_phase(raw)==-1

@pytest.mark.parametrize('field,value',[('coordinate_Wick','I'),('wilson_word',['a1','a2','a3']),
    ('Dyson_external_coefficient',-1),('raw_Euler_count',1),('one_pi_includes_internal_loop_signs',False)])
def test_unsupported_sign_contract_is_rejected(field,value):
    raw=source();raw['physical_convention'][field]=value
    with pytest.raises(ValueError,match='contract_mismatch'):validate_contract(raw,required=True)

def test_missing_native_source_contract_is_rejected():
    with pytest.raises(ValueError,match='requires_physical'):validate_contract({'vertices':[]},required=True)

def test_tree_native_C_has_the_trusted_positive_sign():
    c=compute_c_tgv()
    full=expression(_combine_c_tgv(c))
    # At eps=0 and A_ref=1, the action-derived coefficient is +CF*CA*alphaS^2/(4pi^2).
    full=full.subs({s.Symbol('eps'):0,s.Symbol('I0'):0,s.Symbol('v14'):-s.Rational(1,2)})
    expected=s.Symbol('CA')*s.Symbol('CF')*s.Symbol('alphaS')**2/(4*s.pi**2)
    assert s.simplify(full-expected)==0

def test_FGV_native_C_matches_the_Euclidean_quartic_action():
    c=compute_c_fgv();full=expression(_build_c_fgv_block(c)['combined'])
    full=full.subs({s.Symbol('eps'):0,s.Symbol('I0'):0,s.Symbol('v14'):-s.Rational(1,2)})
    expected=-s.Symbol('CA')**2*s.Symbol('CF')*s.Symbol('alphaS')**3/(16*s.pi**3)
    assert s.simplify(full-expected)==0
    assert coordinate_phase(4,4,1,vertex4=1)=='-1'
