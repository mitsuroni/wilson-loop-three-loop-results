"""Index actual mathematical objects and program roles, not just a final table."""
from prepare import HERE,ROOT,RUNTIME,read,put,sha

R1='integrals/one-loop/';R2='integrals/two-loop/'
M=R2+'sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/'
Q=R2+'sources/workflow/checks/4e04-high-order/runs/20260926T084846954927Z/'
T=R2+'sources/workflow/checks/tgv04-high-order/runs/positive-epsilon/'
def ref(path,object,role,**kwargs):
    assert (RUNTIME/path).is_file(),path
    return {'path':path,'sha256':sha(RUNTIME/path),'object':object,'role':role,**kwargs}

def main():
    for entry in read(HERE/'entries.json'):
        oid=entry['orbit_id'];own=ROOT/entry['path'];index=read(own/'data/index.json');two=entry['loop_directory']=='two-loop'
        primitives=[];codes=[];comparison=[]
        if not two:
            codes.append(ref(R1+'reproduce.py','original_source_checks; derive_coefficients','Rebuild source geometry and actual evaluated Gamma/M generating formulas before reference comparison.'))
            historical=read(own/'input.json')['historical_result'];rel=historical['path'].split('low-loop-runtime/',1)[1]
            comparison.append(ref(rel,'historical complete I','Post-generation holdout comparison, not an integrated subprimitive.'))
            index['external_integrated_primitives_required']=False
        else:
            codes.append(ref(R2+'family_engines.py','reconstruct','Expand declared evaluated mother functions and Gamma/cosh factors. TGV04 instead adopts assembled I coefficients.'))
            codes.append(ref(R2+'reproduce.py','main; exact HLOG comparison','Generate or adopt the stated coefficients, then compare exact representations and previous reports.'))
            codes.append(ref(R2+'sources/workflow/checks/single-opposite/hlog_algebra.py','exact rational shuffle algebra','Function-coordinate arithmetic; not a new parameter integration.'))
            needed=[]
            if oid in ('4e-02','4e-05','3e-3g-03'):needed=['01','11','21','03','31','13','41','23','05']
            elif oid=='4e-03':needed=['01']
            elif oid=='4e-06':needed=['01','11','21','03']
            elif oid=='4e-04':needed=['01','11','21','03','02','12','13','22','23','04','05']
            elif oid=='3e-3g-04':needed=['01','11','21','13','31','12','23','41','22']
            elif oid=='2e-Pi1-02':needed=['01','11','21','03','31','13']
            for k,j in needed:
                primitives.append(ref(M+f'M_{k}_{j}.txt',f'M_{{{k},{j}}} exact HLOG coordinates','Adopted evaluated moment dictionary; the kernel, branch and finite moment budget are stated in the derivation.'))
            if oid=='4e-04':
                for n in range(3):primitives.append(ref(Q+f'Q{n}_plus.txt',f'bare crossing-even Q_{n}','Already integrated two-variable coefficient, without physical prefactor; independent Q-to-I assembly control.',epsilon_order=n))
                codes.append(ref(R2+'sources/workflow/checks/4e04-high-order/compute.mpl','historical HyperInt Q integration','Original integrator and fibration conversion; old absolute CAS setup retained, not invoked by the Python adapter.'))
                index['N_to_M_relation']='N_(k,j)=M_(k,j+1)-ell*M_(k,j); not an extra independent integration file.'
            if oid=='3e-3g-04':
                for script,output,label in [('easy1','D2','easy jet D2'),('easy2','D3','easy jet D3'),('hard1','H1','hard jet H1'),('hard2','H2','hard jet H2'),('recoil1','O1','Euler-subtracted recoil jet O1'),('recoil2','O2','Euler-subtracted recoil jet O2')]:
                    primitives.append(ref(T+output+'.txt',label,'Historical evaluated subintegral. Complete kernels, integration order and exact HLOG value are exposed in the article and appendices.'))
                    codes.append(ref(T+script+'.mpl',label+' integration','Historical Maple/HyperInt integration program, not currently called by the supported Python replay.'))
                for q in (1,2):
                    primitives.append(ref(T+f'I{q}.txt',f'assembled normalized I_{q}','Full-I historical HLOG output for representation comparison; not a D/H/O subintegral.'))
                    codes.append(ref(T+f'assemble{q}.mpl',f'D/H/O to I_{q}','Historical prefactor convolution and assembly. The new Python adapter does not rerun this convolution.'))
                primitives.append(ref(R2+'sources/workflow/checks/tgv04-high-order/compact.json','already assembled normalized I1/I2','Adopted full-I input used by the actual Python replay; not a new integration or live D/H/O assembly.'))
                primitives.append(ref(R2+'sources/benchmarks/2ne/families/single-opposite/basis-functions.json','functions.tauDq4.hlog_terms','Adopted exact HLOG coordinates for tau_D4.'))
            group='pure' if oid.startswith('4e-') else 'self-energy' if oid.startswith('2e-') else 'tgv'
            comparison.append(ref(R2+f'sources/workflow/checks/two-loop-report-20261015/{group}-data.json','previous family coefficient report','Post-generation or post-adoption comparison only.'))
            source=read(own/'input.json')['family_catalog_entry']['canonical_source']
            comparison.append(ref(R2+source['path'],'canonical complete I','Accepted whole-I target, separate from primitive data.'))
        index.update({'path_base':'All paths in adopted_primitives, assembly_code and comparison_targets are relative to runtime_base.','adopted_primitives':primitives,'assembly_code':codes,'comparison_targets':comparison,
            'scope':'Object roles are separate: original initial roots are bound by ../input.json, adopted subintegrals and whole-I inputs are labelled here, comparison targets are not generation inputs. A batch may load additional moment files needed by its other members.'})
        put(own/'data/index.json',index)
    print({'object_level_data_indices':14,'frozen_objects_modified':0})
if __name__=='__main__':main()
