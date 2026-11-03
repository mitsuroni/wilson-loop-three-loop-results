"""Verify the single manuscript, diagram table, and preserved evidence."""
from pathlib import Path
import hashlib
import json
import re
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PI1=ROOT/'benchmarks/self-energy/pi1'
DOC=ROOT/'workflow/checks/self-energy-pi1-baselines/history/benchmarks/self-energy/pi1/theory.md'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    text=DOC.read_text(encoding='utf-8')
    assert '\u200b' not in text
    assert [int(n) for n in re.findall(r'\\tag\{(\d+)\}',text)]==list(range(1,27))
    for unwanted in ('I_g','I_k',r'I_{\rm full}',r'I_{\rm paper}',r'I_{\partial}'):
        assert unwanted not in text,unwanted
    assert 'theory-r001.md' not in text and 'theory-r002.md' not in text
    assert r'(1+2\varepsilon)df+(1-2\varepsilon)uw' in text
    assert r'16\varepsilon(1-2\varepsilon)' in text
    assert r'\varepsilon^3(1-2\varepsilon)(1+\varepsilon)' in text
    checks=[]
    e,z,A,alpha=s.symbols('eps z1 A alpha_s',positive=True)
    CA,CF,TF,nf=s.symbols('CA CF TF nf')
    symbols={'eps':e,'z1':z,'A':A,'alpha_s':alpha,'CA':CA,'CF':CF,'TF':TF,'nf':nf,'Pi':s.pi,'gamma':s.gamma}
    def parse(value):
        return s.sympify(value.replace('Gamma[1-eps]','gamma(1-eps)').replace('(-2*v14+I0)','A').replace('^','**'),locals=symbols)
    def equal(name,value):
        residual=s.factor(s.cancel(s.powsimp(s.expand_power_base(value,force=True),force=True)))
        assert residual==0,(name,residual)
        checks.append({'name':name,'passed':True})
    B=((10-6*e)*CA-8*(1-e)*TF*nf)/(3-2*e)
    rows={}
    for line in text.splitlines():
        if re.match(r'\| [24]e-Pi1-\d\d \|',line):
            parts=[p.strip() for p in line.split('|')[1:-1]]
            rows[parts[0]]=parts[1:]
    assert len(rows)==13
    hashes=[]
    for loops,family in ((2,'2e-Pi1'),(3,'4e-Pi1')):
        source=PI1/f'sources/{family}-r001.json'
        data=json.loads(source.read_text(encoding='utf-8'))
        hashes.append({'path':str(source.relative_to(ROOT)),'sha256':sha(source)})
        for orbit in data['orbits']:
            name=orbit['orbit_id'];j=orbit['J'];geo=j['geometry_polynomials']
            expected_row=[', '.join(geo[0]['s_vars']),', '.join(geo[1]['s_vars']) if len(geo)>1 else '—',
                          ' → '.join(j['iterated_order_outer_to_inner'])]
            assert rows[name]==expected_row,name
            checks.append({'name':name+' endpoints and ordered domain table','passed':True})
            f,d,u,w=[parse(geo[0][key]) for key in ('f','d','u','w')]
            N=(1+2*e)*d*f+(1-2*e)*u*w
            equal(name+' complete kernel numerator',N-parse(j['se_numerator']))
            expectedJ=N/(2*(1+e))*f**(-2+2*e)
            for bare in geo[1:]:expectedJ*=parse(bare['d'])*parse(bare['f'])**(-1+e)
            equal(name+' single complete J',parse(j['integrand'])-expectedJ)
            TW=parse(orbit['color_contraction']['wilson_trace'])
            C=(-1)**(loops-2)*TW*B/(2**(loops+2)*e*(1-2*e))*(alpha/s.pi)**loops*(s.pi**e*s.gamma(1-e))**loops*A**(loops*e)
            equal(name+' common C',parse(orbit['C']['combined'])/C-1)
            if name=='2e-Pi1-01':
                I=C*z**(2*e)/(4*e**2*(1+e))
                final=CF*B/64*(alpha/s.pi)**2*s.gamma(1-e)**2*(s.pi*A*z)**(2*e)/(e**3*(1-2*e)*(1+e))
                equal('Eq1 to Eq22 unique I normalization',I/final-1)
    migration=json.loads((HERE/'navigation.json').read_text(encoding='utf-8'))
    preserved=[]
    for record in migration['historical_files']:
        path=ROOT/record['path']
        if record['path']=='benchmarks/self-energy/pi1/README.md':path=ROOT/migration['readme_snapshot']
        assert sha(path)==record['sha256'],record['path']
        preserved.append(record['path'])
    original=json.loads((ROOT/'workflow/checks/self-energy-pi1/publication.json').read_text(encoding='utf-8'))
    assert original['files']['benchmarks/self-energy/pi1/README.md']==sha(ROOT/migration['readme_snapshot'])
    readme=(PI1/'README.md').read_text(encoding='utf-8')
    assert '(theory.md)' in readme
    assert not any(name in readme for name in ('theory-r001.md','theory-r002.md','paper-check-r001.md','legacy-comparison-r001.md'))
    links=0
    for path in (DOC,PI1/'README.md',PI1/'history.md',HERE/'README.md'):
        body=path.read_text(encoding='utf-8')
        for match in re.finditer(r'\]\(([^)]+)\)',body):
            link=match.group(1)
            if link.startswith(('https://','http://','#')):continue
            target=((PI1 if path==DOC else path.parent)/link).resolve()
            # checks.json is this checker's output, committed only after success.
            if target!=(HERE/'checks.json').resolve():
                assert target.exists(),(path,link)
            links+=1
    report={'canonical_document':str(DOC),'sha256':sha(DOC),'numbered_equations':26,
            'one_I_definition':True,'separate_I_components':False,'orbits':13,
            'mathematical_connections':len(checks),'checks':checks,'preserved_historical_files':preserved,
            'input_hashes':hashes,'local_links_verified':links,'production_math_changed':False,
            'old_publication_manifest_changed':False,'formal':'not_run',
            'limitation':'Document consolidation and exact formula connections; no new integration, production or formal proof claim.'}
    (HERE/'checks-replay.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f'{len(checks)} formula/table connections; 26 equations; one I; {len(preserved)} historical artifacts preserved')

if __name__=='__main__':main()
