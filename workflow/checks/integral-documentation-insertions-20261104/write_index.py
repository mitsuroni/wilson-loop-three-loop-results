"""Append the15-entry reader navigation without replacing earlier39 articles."""
from prepare import HERE,ROOT,read
def main():
    file=ROOT/'integrals/three-loop/README.md';text=file.read_text(encoding='utf-8-sig');marker='\n## Two-point and exchange insertions\n'
    if marker in text:text=text.split(marker)[0].rstrip()+'\n'
    lines=[marker,'This next batch contains15 sources with complete scalar/geometric Laurent coefficients and external colour contributions. Read the [shared insertion mathematics](../families/three-loop-insertions/README.md) and [scientific audit / execution scope](../../workflow/checks/integral-documentation-insertions-20261104/README.md).','','| Source | Route | Mathematics and reproduction |','|---|---|---|']
    for e in read(HERE/'entries.json'):
        oid=e['orbit_id'];lines.append(f"| [{oid}]({oid}/README.md) | {e['family']} | [Derivation]({oid}/derivation.md), [reproduce]({oid}/reproduce.md), [validation]({oid}/validation.md) |")
    lines.extend(['','Every adapter runs the same15-entry batch. Source-defined signs, the full transverse tensor, endpoint contributions and all external operators are retained. Pi2 geometric I is not itself a complete colour channel. The complete full-D internal Pi2 operator and evaluated moment functions are explicit adopted inputs.',''])
    file.write_text(text+'\n'.join(lines),encoding='utf-8');print({'new_entries':15,'earlier39_entries_retained':True})
if __name__=='__main__':main()
