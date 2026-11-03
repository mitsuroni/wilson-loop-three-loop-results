"""Check the integrated full-epsilon Z sum against a separate 1D quadrature."""
import json
import mpmath as mp
from formulas import run,numerical_full
def main():
    mp.mp.dps=50
    ep=mp.mpf("0.3");N=32
    expected=numerical_full(ep,1,50)["Z"]
    B=mp.beta(2-ep,1-ep);total=mp.mpf(0)
    for k in range(1,N+1):
        F=B/(k+1+ep)*mp.hyper([1,1,2-ep],[k+2+ep,3-2*ep],1)
        total-=mp.rf(-ep,k)/(mp.factorial(k)*(k+2*ep))*F
    # 1/(1-x+x*r)<=1/r; sum_{k>N} -(-e)_k/k! = (1-e)_N/N!.
    tail=B*mp.rf(1-ep,N)/mp.factorial(N)/((N+1+2*ep)*(N+1+ep))
    residual=expected-total
    assert 0<residual<tail
    out={"passed":True,"eps":str(ep),"terms":N,"partial_sum":str(total),
        "one_dimensional_quadrature":str(expected),"difference":str(residual),
        "analytic_positive_tail_upper_bound":str(tail),
        "bound_note":"analytic bound for real 0<eps<1; displayed numeric evaluation is not interval arithmetic",
        "formal":"not_run"}
    (run()/"series-check.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))
if __name__=="__main__":main()
