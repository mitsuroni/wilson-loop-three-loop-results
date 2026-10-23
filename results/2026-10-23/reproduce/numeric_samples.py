"""Evaluate delivered moment formulas by their real, removable-kernel integrals.

These are reproducible representation samples, not independent QCD integrals.
"""
from pathlib import Path
import argparse, hashlib, json
import mpmath as mp
import sympy as s

def moment(i, j, z):
    rho = mp.log(z)
    r = z/(1+z)
    def kernel(t):
        if t == 0 or t == 1:
            return mp.mpf('0')  # isolated integration endpoint, not its limit
        w = mp.log(t)+mp.log1p(-t)
        if t == r:
            return -w**i/(r*(1-r)) if j == 1 else mp.mpf('0')
        ell = rho+mp.log1p(-t)-mp.log(t)
        return w**i*ell**j/(t-r)
    return mp.quad(kernel, [0, r/2, r, (1+r)/2, 1])

def evaluate(expr, z, scale_L):
    raw = s.sympify(expr)
    vals = {'rho':mp.log(z), 'L':mp.mpf(scale_L)}
    for symbol in raw.free_symbols:
        name = str(symbol)
        if name.startswith('M_'):
            _, i, j = name.split('_'); vals[name] = moment(int(i), int(j), z)
        elif name.startswith('N_'):
            _, i, j = name.split('_'); i, j = int(i), int(j)
            vals[name] = moment(i, j+1, z)-mp.log(z)*moment(i, j, z)
    syms = sorted(raw.free_symbols, key=str)
    return s.lambdify(syms, raw, modules='mpmath')(*[vals[str(x)] for x in syms])

def main():
    p = argparse.ArgumentParser(); p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args = p.parse_args(); out = args.output_dir.resolve()
    source = out/'coefficients.json'; data=json.loads(source.read_text(encoding='utf-8'))
    rows=[]; checks=[]
    for color, rec in data['channels'].items():
        for ztext in ['0.5','1','2']:
            for scale_L in ['0','1']:
                values=[]
                for dps in [40,60]:
                    with mp.workdps(dps): values.append(evaluate(rec['F'],mp.mpf(ztext),scale_L))
                with mp.workdps(60):
                    delta=abs(values[0]-values[1]); assert delta < mp.mpf('1e-32'), (color,ztext,scale_L,delta)
                    row={'channel':color, 'z':ztext,'L':scale_L,'F':mp.nstr(values[1],40),'precision_change_absolute':mp.nstr(delta,5)}
                    rows.append(row); checks.append({'name':f'{color} z={ztext} L={scale_L} 40/60 digit stability','passed':True})
    for color in data['channels']:
        for scale_L in ['0','1']:
            pair=[x for x in rows if x['channel']==color and x['L']==scale_L and x['z'] in ['0.5','2']]
            with mp.workdps(60): delta=abs(mp.mpf(pair[0]['F'])-mp.mpf(pair[1]['F']))
            assert delta < mp.mpf('1e-35')
            checks.append({'name':f'{color} L={scale_L} reciprocal-z sample','passed':True})
    result={'status':'passed','date':'2026-10-23','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'scope':'Direct quadrature of the defined M/N moments in F; representation evaluation only, not independent integration of the original diagrams.',
            'mpmath_version':mp.__version__,'working_precisions':[40,60],'reported_digits':40,'samples':rows,'checks':checks}
    (out/'numeric-samples.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'samples':len(rows),'checks':len(checks)}))

if __name__=='__main__': main()
