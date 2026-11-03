"""Reuse the checked general reflection map on the three T2[1] ordered inputs."""
from pathlib import Path
import json,sys,datetime,importlib.util,copy
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("t2empty_bridge",ROOT/"workflow/checks/tgv-t2-empty/check_math.py")
bridge=importlib.util.module_from_spec(spec);spec.loader.exec_module(bridge)
read,dump,sha=bridge.read,bridge.dump,bridge.sha
PAIR={"16":"32","36":"12","37":"11"}
def run():return ROOT/read(HERE/"current-run.json")["run"]
def compute():
    raw=bridge.records(bridge.SOURCE);roots=bridge.records(bridge.PREPARED)
    catalog=read(ROOT/"benchmarks/tgv/families/catalog-r001.json")
    members={oid for g in catalog["groups"] if g["group_id"] in ("t2-opposite-13","t2-opposite-24") for oid in g["mixed_members"]}
    assert members=={"5e-3g-"+x for x in PAIR}
    manifest=read(bridge.PREPARED.parent/"manifest.json");assert sha(bridge.PREPARED)==manifest["sha256"]
    # Find matches among all three source roots rather than accept labels alone.
    candidates={};rows=[]
    for target in PAIR:
        found=[]
        for source in ("11","12","32"):
            try:r=bridge.check_pair(raw[target],raw[source],roots[target],roots[source])
            except AssertionError:continue
            found.append(source)
        assert found==[PAIR[target]],(target,found)
        candidates[target]=found
        rows.append(bridge.check_pair(raw[target],raw[PAIR[target]],roots[target],roots[PAIR[target]]))
    rejected=[]
    for name in ("missing_alpha_permutation","missing_order_reversal","wrong_P_sign","wrong_C_factor","wrong_pair"):
        t=copy.deepcopy(raw["37"]);rt=copy.deepcopy(roots["37"]);s=raw["11"];rs=roots["11"];kwargs={}
        if name=="missing_alpha_permutation":kwargs["swap_alpha"]=False
        elif name=="missing_order_reversal":kwargs["reverse_order"]=False
        elif name=="wrong_P_sign":
            for d in (t,rt):
                old=d["J_raw_s"]["tgv"]["pnum"];d["J_raw_s"]["tgv"]["pnum"]="-("+old+")"
        elif name=="wrong_C_factor":rt["C"]["expression"]="2*("+rt["C"]["expression"]+")"
        elif name=="wrong_pair":s=raw["12"];rs=roots["12"]
        try:bridge.check_pair(t,s,rt,rs,**kwargs)
        except AssertionError:rejected.append(name)
        else:raise AssertionError("bad bridge accepted: "+name)
    z=bridge.sp.Symbol("z1")
    gram=bridge.sp.Matrix([[0,z,-1-z,1],[z,0,1,-1-z],[-1-z,1,0,z],[1,-1-z,z,0]])
    assert gram==gram.extract([1,0,3,2],[1,0,3,2])
    prep=read(bridge.PREPARED)
    assert prep["I_definition"]=="I=sym/32*Exp[-3*EulerGamma*eps]*Gamma[1-eps]*Gamma[2-2*eps]*(J(z,eps)+z^(3*eps)*J(1/z,eps))/(1+z^(3*eps))"
    return {"passed":True,"pairs":rows,"unique_candidate_matches":candidates,"negative_cases_rejected":rejected,
      "Gram_and_reference_scale_unchanged":True,"z_transformation":"z -> z",
      "full_epsilon_relations":["I16=I32","I36=I12","I37=I11","I36+I37=I16"],
      "domain":"z>0; first 1/2<Re(eps)<1, then meromorphic epsilon continuation","formal":"not_run"}
def main():
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    r=HERE/"runs"/stamp;r.mkdir(parents=True);dump(HERE/"current-run.json",{"run":r.relative_to(ROOT).as_posix()})
    result=compute();dump(r/"reflection.json",result)
    prior=ROOT/read(ROOT/"workflow/checks/tgv-t1-opposite/current-delivery.json")["delivery"]
    dump(r/"source-delivery.json",{"delivery":prior.relative_to(ROOT).as_posix(),"manifest_sha256":sha(prior/"manifest.json"),
       "evidence_run":read(prior/"manifest.json")["evidence_run"]})
    dump(r/"environment.json",{"python":sys.version,"sympy":bridge.sp.__version__,"source_sha256":sha(bridge.SOURCE),"prepared_sha256":sha(bridge.PREPARED)})
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
