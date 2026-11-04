"""New Pi2 operator/colored remainder checks using bound, already solved geometry."""
from pathlib import Path
import importlib.util,json,re,hashlib
import mpmath as mp
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi2'
geometry_path=ROOT/'workflow/checks/self-energy-pi1-chain-solve/numeric.py'
spec=importlib.util.spec_from_file_location('chain_numeric_reuse',geometry_path)
geo=importlib.util.module_from_spec(spec);spec.loader.exec_module(geo)
mp.mp.dps=80
records=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
    package=load(BASE/'prepared/2e-Pi2.json');expansion=load(BASE/'operator-expansion-r001.json')
    ep=s.Symbol('eps');U=s.gamma(1+ep)**2*s.gamma(1-ep)*s.gamma(1-3*ep)/(s.gamma(1-2*ep)**2*s.gamma(1+2*ep))
    phi={}
    for name,c in package['orbits'][0]['operator']['channels'].items():
        A=s.sympify(c['A'].replace('^','**'),locals={'eps':ep});B=s.sympify(c['B'].replace('^','**'),locals={'eps':ep})
        phi[name]=s.lambdify(ep,4*ep**2*(A+B/U),'mpmath')
    allvalues={}
    for zz in ['.5','1','2']:
        allvalues[zz]=geo.values(mp.mpf(zz))
        for name in ['M21','M03']:
            records.append({'type':'existing defining M moment','z':zz,'name':name,'value':str(allvalues[zz][name]),'passed':True})
    print('6 defining moment evaluations completed.',flush=True)
    denominators=[128,256,512,1024,2048]
    for name in ['AA','Af','Ff']:
        cs={int(k):s.lambdify((),s.sympify(v['expression']),'mpmath')() for k,v in expansion['channels'][name].items()}
        errors=[]
        for den in denominators:
            e=mp.mpf(1)/den;errors.append(abs(phi[name](e)-mp.fsum(v*e**k for k,v in cs.items())))
        ratio=errors[-1]/errors[-2]
        assert abs(ratio-mp.mpf(1)/32)<mp.mpf('.003'),(name,ratio)
        records.append({'type':'phi expansion O(eps^5)','channel':name,'denominators':denominators,'errors':list(map(str,errors)),'last_ratio':str(ratio),'passed':True})
    for key in ['01','02']:
        result=load(BASE/f'results/2e-Pi2-{key}/r001.json')
        for zz in ['.5','1','2']:
            z=mp.mpf(zz);v=allvalues[zz];args=[v[n] for n in ['L','M01','M11','M21','M03']]
            for channel in result['color_channels']:
                name=channel['channel'];fs=geo.coeff_functions(channel['laurent']['coefficients']);cs={k:f(*args) for k,f in fs.items()}
                errors=[]
                for den in denominators:
                    e=mp.mpf(1)/den
                    errors.append(abs(phi[name](e)*geo.I(key,e,z)-mp.fsum(v*e**k for k,v in cs.items())))
                ratios=[errors[k+1]/errors[k] for k in range(len(errors)-1)]
                record={'type':'full color I O(eps) remainder','orbit':key,'channel':name,'z':zz,'denominators':denominators,'errors':list(map(str,errors)),'ratios':list(map(str,ratios)),'passed':bool(abs(ratios[-1]-mp.mpf('.5'))<mp.mpf('.03'))}
                if not record['passed']:
                    (HERE/'failed-remainder-diagnostic.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
                assert record['passed'],(key,name,zz,ratios)
                records.append(record)
            print('color remainders '+key+' z='+zz+' passed',flush=True)
    # Reconstruct the actual full-color prepared C/J, including Euler only once.
    for o in package['orbits']:
        text=o['C']['combined'].replace('(-2*v14+I0)','scale')
        text=re.sub(r'Gamma\[([^\]]+)\]',r'gamma(\1)',text).replace('^','**')
        ep,ca,cf,tf,nf,alpha,scale=s.symbols('eps CA CF TF nf alpha_s scale')
        obj=s.sympify(text,locals={'eps':ep,'CA':ca,'CF':cf,'TF':tf,'nf':nf,'alpha_s':alpha,'scale':scale,'Pi':s.pi,'EulerGamma':s.EulerGamma,'gamma':s.gamma,'exp':s.exp})
        fn=s.lambdify((ep,ca,cf,tf,nf,alpha,scale),obj,'mpmath')
        e=mp.mpf('.13');z=mp.mpf('.7');al=mp.mpf('.2');A=mp.mpf('2.6');key=o['orbit_id'][-2:]
        for flavors in [0,4,6]:
            C=fn(e,3,mp.mpf(4)/3,mp.mpf('.5'),flavors,al,A)
            obtained=mp.mpf(o['I_target']['sym'])/2*C*(geo.J(key,e,z)+z**(3*e)*geo.J(key,e,1/z))/((al/mp.pi)**3*(mp.pi*A)**(3*e)*(1+z**(3*e)))
            colors={'AA':mp.mpf(12),'Af':mp.mpf(2)*flavors,'Ff':mp.mpf(8)/9*flavors}
            expected=mp.fsum(colors[c]*phi[c](e)*geo.I(key,e,z) for c in ['AA','Af','Ff'])
            error=abs(obtained/expected-1);assert error<mp.mpf('1e-60')
            records.append({'type':'actual prepared C/J full-color assembly','orbit':key,'nf':flavors,'relative_error':str(error),'passed':True})
    files=[geometry_path,BASE/'prepared/2e-Pi2.json',BASE/'operator-expansion-r001.json',*[BASE/f'results/2e-Pi2-{key}/r001.json' for key in ['01','02']]]
    report={'count':len(records),'checks':records,'mpmath_dps':80,'dependencies':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in files],
        'scope':'New Pi2 scalar series, colored Laurent remainders and prepared full-color reconstruction; reused geometry is bound by exact kernel/domain checks.',
        'new_original_kernel_quadrature':False,'formal':'not_run'}
    (HERE/'numeric-checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'{len(records)} numerical controls/moment evaluations passed.',flush=True)
if __name__=='__main__':main()
