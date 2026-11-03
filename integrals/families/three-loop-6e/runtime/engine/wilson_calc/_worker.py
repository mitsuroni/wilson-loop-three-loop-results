"""Internal subprocess protocol; accepts only closed scalar AST payloads."""
import json
import sys
import sympy as sp
from .algebra import ast_to_sympy
from .assembly import expand_scalar


def main():
    from interfaces.validate import load
    value=load(sys.argv[1])
    eps=sp.Symbol('eps');z=sp.Symbol('z1',positive=True)
    j=ast_to_sympy(value['J'],{'eps':eps,'z1':z})
    coefficients,lower=expand_scalar(j,value['normalization'],value['quantity'],value['through'],eps,z,value['loop_order'])
    print(json.dumps({'format':'wilson-scalar-worker/0.1','coefficients':coefficients,'lower_bound':lower}))


if __name__=='__main__':main()
