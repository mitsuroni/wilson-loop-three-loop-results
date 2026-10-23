"""Apply only the independently derived path-order bridge to the final sum."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'physical-bridge.json').read_text(encoding='utf-8'))
assert d['passed'] and not d['RG_used_to_derive_sign']
assert d['weights']=={'Pi1':1,'chain':1,'Pi2':1,'legs':-1,'Gamma':-1}
source=(HERE/'assemble.mpl').read_text(encoding='utf-8')
old='Raw[q]:=FB(add(parts[p,q],p in partnames)):'
new='Raw[q]:=FB(parts["Pi1",q]+parts["chain",q]+parts["Pi2",q]-parts["legs",q]-parts["Gamma",q]):'
assert old in source
source=source.replace(old,new).replace('assembly-coefficients.txt','physical-coefficients.txt').replace('"assembly.m"','"physical.m"').replace('CHANNEL_ASSEMBLY_COMPLETED','PHYSICAL_CHANNEL_COMPLETED')
source=source.replace('Post[q]-add(parts["Gamma",j]*Scale(3,q-j),j=-6..q)','Post[q]+add(parts["Gamma",j]*Scale(3,q-j),j=-6..q)')
(HERE/'physical.mpl').write_text('# Physical weights are derived in physical-bridge.json, independently of RG.\n'+source,encoding='utf-8')
print('Physical assembly script generated from independently derived bridge.')
