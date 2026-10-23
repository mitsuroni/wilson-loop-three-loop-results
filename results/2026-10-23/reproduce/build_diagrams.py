"""Export the 31 selected matter diagrams and their 35 channel contributions.

Reads scientific inputs and a completed rg-replay; writes only diagram-ledger.json
and diagrams.md in --output-dir. Reuses the established parser for exact algebra,
without integration, selector updates, or a new mixed-argument Hlog sum bridge.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
POWERS = list(range(-6, 1))
COLORS = {'Af': 'CF*CA*TF*nf', 'Ff': 'CF^2*TF*nf', 'ff': 'CF*(TF*nf)^2'}
CHECKER = 'workflow/checks/matter-rg-20261017/check.py'
SOURCES: dict[str, dict] = {}
CACHE: dict[str, object] = {}
ROOT = Path()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def field_hash(value):
    return sha(json.dumps(value, sort_keys=True, ensure_ascii=False,
                          separators=(',', ':')).encode('utf-8'))


def source(path, expected=None, strict=False, role='definition or dependency'):
    path = str(path).replace('\\', '/')
    full = (ROOT / path).resolve()
    if not full.is_relative_to(ROOT):
        raise ValueError(f'Source outside the read-only root: {path}')
    data = full.read_bytes()
    actual = sha(data)
    entry = SOURCES.setdefault(path, {'path': path, 'sha256': actual,
                                     'bytes': len(data), 'roles': []})
    assert entry['sha256'] == actual, ('concurrent source change', path)
    if role not in entry['roles']:
        entry['roles'].append(role)
    result = {'path': path, 'sha256': actual}
    if expected:
        result.update(declared_sha256=expected, declared_hash_matches=(actual == expected))
        if strict:
            assert actual == expected, ('source binding mismatch', path, expected, actual)
    return result


def read(path, expected=None, strict=False, role='definition or dependency'):
    ref = source(path, expected, strict, role)
    if path not in CACHE:
        CACHE[path] = json.loads((ROOT / path).read_text(encoding='utf-8-sig'))
    return CACHE[path], ref


def ptr(data, pointer):
    for token in pointer.lstrip('/').split('/') if pointer else []:
        token = token.replace('~1', '/').replace('~0', '~')
        data = data[int(token)] if isinstance(data, list) else data[token]
    return data


def field_ref(path, pointer, value):
    return {**source(path, role='selected result or exact representation'),
            'json_pointer': pointer, 'field_sha256': field_hash(value)}


def all_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for x in value.values():
            yield from all_strings(x)
    elif isinstance(value, list):
        for x in value:
            yield from all_strings(x)


def metadata(data):
    """Source-scoped definitions, kept as recorded; not a compact-basis rewrite."""
    exact = {'normalization', 'domain', 'initial_region', 'definition', 'functions',
             'coefficient_notation', 'variables', 'Hlog_variables', 'Hlog_convention',
             'named_factors', 'coefficient_field', 'expression_language',
             'Hlog_expression_language', 'Hlog_language', 'regularization',
             'reference_regularization', 'pole_budget', 'highest_pole_order',
             'vanishing_higher_poles', 'zero_powers', 'range', 'laurent_remainder',
             'z_dependence', 'physical_convention', 'sign_convention', 'formal',
             'limitations', 'completion_semantics', 'native_verification_scope',
             'representation_revision', 'representation_verification',
             'color_coefficient_semantics', 'color_factors', 'R', 'I', 'I_geom',
             'J', 'all_epsilon', 'I_expression_language', 'I_trigonometric'}
    return {k: v for k, v in data.items() if k in exact or k.endswith('_definition')}


def resolve_related(path, base):
    if (ROOT / path).is_file():
        return path
    candidate = (ROOT / base / path).resolve()
    if candidate.is_relative_to(ROOT) and candidate.is_file():
        return candidate.relative_to(ROOT).as_posix()
    return None


def source_contracts(data, oid, base):
    out = []
    inp = data.get('input', {})
    for key, hash_key in [('package', 'package_sha256'),
                          ('raw_package', 'raw_package_sha256'),
                          ('native_package', 'native_package_sha256')]:
        if key not in inp:
            continue
        path = resolve_related(inp[key], base)
        assert path, ('missing native contract', inp[key])
        package, ref = read(path, inp.get(hash_key), True, 'native raw/prepared C/J contract')
        idx = next(i for i, x in enumerate(package['orbits']) if x['orbit_id'] == oid)
        orbit = package['orbits'][idx]
        out.append({**ref, 'json_pointer': f'/orbits/{idx}',
                    'kind': key, 'geometry': package.get('geometry'),
                    'normalization_contract': package.get('normalization_contract'),
                    'fields': {k: {'json_pointer': f'/orbits/{idx}/{k}',
                                   'field_sha256': field_hash(orbit[k])}
                               for k in ['C', 'J', 'J_raw_s', 'J_ibp', 'J_sector',
                                         'branch_prescription', 'normalization',
                                         'operator', 'color_contraction', 'I_target']
                               if k in orbit},
                    'raw_definitions': {k: orbit[k] for k in
                                        ['branch_prescription', 'normalization',
                                         'operator', 'color_contraction', 'I_target']
                                        if k in orbit}})
    return out


def dependencies(data, base):
    out = []
    for item in data.get('dependencies', []):
        path = resolve_related(item['path'], base)
        if path:
            out.append(source(path, item.get('sha256'), role='source symbol definition dependency'))
    for key in ['derivation', 'representation_verification', 'native_verification', 'evidence']:
        value = data.get(key)
        if isinstance(value, str):
            path = resolve_related(value, base)
            if path:
                out.append(source(path, role=key))
    return out


def coefficient_map(data, path, pointer, parser):
    raw = ptr(data, pointer)
    output = {}
    for q in POWERS:
        key = str(q)
        if key in raw:
            value = raw[key]
            parsed = parser.parse(value).subs(parser.L, parser.rho)
            output[key] = {'raw': value, 'origin': 'explicit source coefficient',
                           'source': field_ref(path, pointer + '/' + key, value),
                           'normalized_expression': str(parsed),
                           'is_zero': parsed == 0}
        else:
            higher = data.get('vanishing_higher_poles', {})
            if key in higher and str(higher[key]) == '0':
                zero_ref = '/vanishing_higher_poles/' + key
            elif q in data.get('zero_powers', []):
                zero_ref = '/zero_powers'
            elif 'highest_pole_order' in data and q < -data['highest_pole_order']:
                zero_ref = '/highest_pole_order'
            else:
                raise ValueError(f'Unexplained missing coefficient: {path}{pointer}/{key}')
            output[key] = {'raw': None, 'origin': 'zero established by source pole metadata',
                           'source': field_ref(path, zero_ref, ptr(data, zero_ref)),
                           'normalized_expression': '0', 'is_zero': True}
    return {'source': field_ref(path, pointer, raw), 'coefficients': output,
            'range': [-6, 0], 'remainder': 'O(eps)',
            'raw_field': raw}


def input_field(data, part, channel):
    if part == 'Pi1':
        return '/laurent/coefficients'
    if part == 'legs':
        return '/laurent_Hlog_coefficients' if 'laurent_Hlog_coefficients' in data else '/laurent_coefficients'
    if part == 'Gamma':
        return '/channels/q/I_coefficients'
    for i, c in enumerate(data['color_channels']):
        if part == 'chain' and c['monomial'] == ('CA*CF*TF*nf' if channel == 'Af' else COLORS[channel]):
            return f'/color_channels/{i}/coefficients'
        if part == 'Pi2' and c['channel'] == channel:
            return f'/color_channels/{i}/laurent/coefficients'
    raise ValueError((part, channel))


def compact_alternatives(data, path, selector, base, part):
    out = []
    if part == 'legs':
        for field in ['laurent_coefficients', 'M_representation', 'second_pole_in_existing_M_basis']:
            if field in data and (field != 'laurent_coefficients' or 'laurent_Hlog_coefficients' in data):
                out.append({'label': field, 'source': field_ref(path, '/' + field, data[field]),
                            'raw': data[field], 'definitions': metadata(data),
                            'validation': data.get('representation_verification'),
                            'status': 'existing source representation; no new conversion performed'})
        for p in selector.get('equivalent_representations', []):
            apath = base + '/' + p['result']
            alt, ref = read(apath, p['sha256'], True, 'selector-listed compact alternative')
            out.append({'label': p.get('scope', alt.get('status')), 'source': ref,
                        'raw': alt, 'status': 'selector-listed equivalent representation',
                        'validation': alt.get('evidence')})
            dependencies(alt, base)
    if part == 'Gamma' and 'preferred_representation' in selector:
        apath = base + '/' + selector['result']
        alt, ref = read(apath, selector['sha256'], True, 'default representation of preferred z result')
        out.append({'label': 'Default r representation of the same quark channel',
                    'source': field_ref(apath, '/channels/q/I_coefficients', alt['channels']['q']['I_coefficients']),
                    'raw': alt['channels']['q']['I_coefficients'], 'definitions': metadata(alt),
                    'validation': data.get('conversion'),
                    'status': 'existing source conversion; no new r-to-z conversion performed'})
    return out


def main():
    global ROOT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    ROOT = args.root.resolve()
    output = args.output_dir.resolve()
    assert output.is_dir(), 'Prepare output-dir and its completed rg-replay first'
    replay = output / 'rg-replay'
    rp = json.loads((replay / 'inputs.json').read_text(encoding='utf-8'))
    rr = json.loads((replay / 'current-coefficients.json').read_text(encoding='utf-8'))
    rc = json.loads((replay / 'checks.json').read_text(encoding='utf-8'))
    assert rc['status'] == 'passed'
    for path, expected in rp['source_sha256'].items():
        source(path, expected, True, 'fresh RG replay input')
    spec = importlib.util.spec_from_file_location('matter_rg_parser', ROOT / CHECKER)
    parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parser)
    source(CHECKER, role='exact parser and rational convolution implementation')
    old, _ = read('workflow/checks/w3-ca-nf/inputs.json')
    old_by_id = {x['orbit_id']: x for x in old['parts']}
    selectors = {}
    for key, path in [('se', 'benchmarks/self-energy/current-r002.json'),
                      ('legs', 'benchmarks/self-energy/tgv-pi1/current.json'),
                      ('Gamma', 'benchmarks/self-energy/gamma1/current.json')]:
        selectors[key], _ = read(path, role='current selector')
    for path in ['benchmarks/self-energy/conventions.md',
                 'benchmarks/2ne/families/single-opposite/basis.md',
                 'benchmarks/2ne/families/single-opposite/basis-functions.json',
                 'benchmarks/theory/montopar/shared-sign-conventions.md',
                 'benchmarks/renormalization/three-loop/assembly-ledger.md',
                 'benchmarks/renormalization/three-loop/matter-channels.md']:
        if (ROOT / path).is_file():
            source(path, role='normalization, notation, or assembly contract')
    diagrams, checks = [], []
    for record in rp['graph_consumption']:
        oid, part, path = record['orbit_id'], record['part'], record['path']
        data, ref = read(path, record['sha256'], True, 'selected native result')
        if part in ['legs', 'Gamma']:
            base = 'benchmarks/self-energy/' + ('tgv-pi1' if part == 'legs' else 'gamma1')
            selected = selectors[part]['orbits'][oid]
        else:
            base = 'benchmarks/self-energy/' + {'Pi1': 'pi1', 'chain': 'pi1-chain', 'Pi2': 'pi2'}[part]
            selected = next(x for x in selectors['se']['entries'] if x['orbit_id'] == oid)
        af_field = input_field(data, part, 'Af')
        assert field_hash(ptr(data, af_field)) == record['field_sha256']
        native_field = '/laurent/coefficients' if part == 'chain' else '/I_geom_laurent/coefficients' if part == 'Pi2' else af_field
        native = coefficient_map(data, path, native_field, parser)
        native['meaning'] = ('unweighted geometry I_geom; channel operators are already included in the consumed channel fields' if part == 'Pi2'
                             else 'unweighted scalar I; chain color multipliers are already included in the consumed channel fields' if part == 'chain'
                             else 'quark-channel scalar I_q after stripping the external color monomial' if part == 'Gamma'
                             else 'native scalar I before the external matter coefficient Bf(eps)')
        deps = dependencies(data, base)
        derivation = f'{base}/workspaces/{oid}/derivation.md'
        if (ROOT / derivation).is_file():
            deps.append(source(derivation, role='per-diagram derivation'))
        entry = {'orbit_id': oid, 'family': part, 'result_source': ref,
                 'revision': data.get('revision'), 'selector_entry': selected,
                 'native_unweighted': native, 'source_metadata': metadata(data),
                 'source_contracts': source_contracts(data, oid, base),
                 'dependencies': deps, 'compact_source_alternatives': compact_alternatives(data, path, selected, base, part),
                 'contributions': {}, 'formal': data.get('formal', 'not_run')}
        channels = ['Af'] + (['ff'] if part == 'chain' else ['Ff'] if part == 'Pi2' else [])
        for channel in channels:
            field = input_field(data, part, channel)
            cmap = coefficient_map(data, path, field, parser)
            replay_record = next(x for x in rr[COLORS[channel]]['sources'] if x['orbit_id'] == oid)
            assert replay_record['path'] == path and replay_record['sha256'] == ref['sha256']
            assert replay_record['field_sha256'] == cmap['source']['field_sha256']
            vals = {int(q): parser.parse(v['normalized_expression']) for q, v in cmap['coefficients'].items()}
            external = str(parser.BF) if part in ['Pi1', 'legs'] else '1'
            weighted = parser.multiply_bf(vals) if part in ['Pi1', 'legs'] else vals
            included = None
            if part in ['chain', 'Pi2']:
                color_record = ptr(data, field.rsplit('/laurent/coefficients', 1)[0] if part == 'Pi2' else field.rsplit('/coefficients', 1)[0])
                included = {k: v for k, v in color_record.items() if k not in ['coefficients', 'laurent']}
            out_coeffs = {}
            for q in POWERS:
                value = weighted[q]
                terms = [{'input_epsilon_power': j, 'rational_weight': str(parser.BFC[q-j])}
                         for j in POWERS if j <= q] if part in ['Pi1', 'legs'] else [{'input_epsilon_power': q, 'rational_weight': '1'}]
                out_coeffs[str(q)] = {'expression': str(value), 'language': 'SymPy/str',
                                      'latex': parser.s.latex(value), 'convolution_terms': terms,
                                      'is_zero': value == 0}
                if channel == 'Af':
                    historical = record['legacy_physical_sign'] * parser.parse(old_by_id[oid]['coefficients'][str(q)])
                    residual = parser.norm(value - historical)
                    assert residual == 0, ('weighted diagram mismatch', oid, q, str(residual))
                    checks.append({'orbit_id': oid, 'channel': channel, 'epsilon_power': q,
                                   'check': 'current native weighted coefficient equals sign-corrected historical assembly input',
                                   'passed': True})
            entry['contributions'][channel] = {'color_monomial': COLORS[channel],
                'consumed_native_field': cmap, 'external_scalar_multiplier': external,
                'fresh_replay_field_binding': replay_record,
                'native_physical_sign': 1, 'legacy_physical_sign_for_Af': record['legacy_physical_sign'],
                'color_operator_already_included': included,
                'coefficient_stage': 'unit-scale physical web contribution; before S3 and coupling counterterms',
                'coefficients': out_coeffs}
        # Definitions retain the original source namespace; normalized expressions use rho.
        expressions = '\n'.join(c['expression'] for cset in entry['contributions'].values() for c in cset['coefficients'].values())
        hsymbols = sorted(set(re.findall(r'\bH_[zr]_[A-Za-z0-9_]+', expressions)))
        entry['normalized_symbol_definitions'] = {
            'rho': 'log(z); the geometry logarithm. Any L in the raw diagram result means log(z), not the renormalization-scale L.',
            'r': '1/(1+z)', 'z': 'positive real v12/v14; z1 in the raw source',
            'H_symbols': {h: {'argument': h.split('_')[1],
                              'word': [-1 if x == 'm1' else int(x) for x in h.split('_')[2:]],
                              'definition': 'Hlog(argument,word) with kernels dt/(t-a); no r-to-z rewrite applied'} for h in hsymbols},
            'zeta': 'Riemann zeta',
            'M_N_and_other_aliases': 'Definitions and alternate expressions are source-scoped in source_metadata, compact_source_alternatives, and dependencies; no global renaming to a new compact basis is assumed.'}
        diagrams.append(entry)
    assert len(diagrams) == 31 and len({d['orbit_id'] for d in diagrams}) == 31
    counts = collections.Counter(c for d in diagrams for c in d['contributions'])
    assert counts == {'Af': 31, 'ff': 2, 'Ff': 2}
    # Exact conservation of selected raw fields and bindings is checked before publication.
    for path, rec in SOURCES.items():
        assert sha((ROOT / path).read_bytes()) == rec['sha256'], ('concurrent source change', path)
    result = {'format': 'wilson-three-matter-diagram-ledger/1', 'date': '2026-10-23',
              'read_only_source_root': str(ROOT), 'colors': COLORS,
              'summary': {'unique_diagrams': 31, 'channel_occurrences': 35,
                          'occurrences_by_channel': dict(counts),
                          'family_counts': dict(collections.Counter(d['family'] for d in diagrams)),
                          'coefficient_range': [-6, 0], 'weighted_coefficient_count': 245,
                          'native_unweighted_coefficient_count': 217,
                          'coefficient_fields_matched_to_fresh_replay': 35,
                          'native_unweighted_zero_origin_counts': dict(collections.Counter(
                              v['origin'] for d in diagrams for v in d['native_unweighted']['coefficients'].values() if v['is_zero'])),
                          'fresh_per_diagram_comparison_count': len(checks)},
              'assembly_contract': {'unit_scale_sum': 'u3_c = sum_d contribution_d,c',
                  'physical_bare_sum': 'hat_w3_c = S3 * u3_c, S3=exp(3*eps*(L+rho/2))+exp(3*eps*(L-rho/2))',
                  'Bf': '-8*(1-eps)/(3-2*eps)',
                  'Bf_series_through_eps6': {str(q): str(v) for q, v in parser.BFC.items()},
                  'extra_sign_C4_crossing_Euler_factors': 'none; already included in each selected native I',
                  'pure_AA_channel': 'outside this three-channel ledger',
                  'epsilon_budget': {'three_loop': 0, 'two_loop': 1, 'one_loop': 2},
                  'sum_bridge_scope': 'Af current weighted diagram coefficients match the historical physical input. The mixed r-to-z full-sum normalization reuses the preserved HyperInt witness; it is not freshly recomputed here.',
                  'sum_bridge_sources': ['workflow/checks/w3-ca-nf/physical.mpl', 'workflow/checks/w3-ca-nf/physical-coefficients.txt', 'workflow/checks/w3-ca-nf/physical-bridge.json'],
                  'formal': 'not_run'},
              'raw_definition_policy': 'Raw values and source-scoped definitions are retained unchanged. Legacy source names and syntax are archival data; they are not a newly standardized compact representation.',
              'parser_normalization': {'source': source(CHECKER),
                  'operations': ['exact rational arithmetic and Bf convolution', 'source hlog_terms use z and zeta powers [2,3,5]',
                                 'Hlog arguments and words preserved as H_z_* or H_r_* symbols',
                                 'pure zero z words become rho^n/n!', 'raw diagram L and log(z) become rho'],
                  'not_performed': ['new integration', 'new mixed-argument Hlog canonicalization', 'new compact-function conversion', 'formal verification']},
              'rg_replay': {f: {'path': 'rg-replay/' + f, 'sha256': sha((replay / f).read_bytes())}
                            for f in ['inputs.json', 'current-coefficients.json', 'checks.json']},
              'diagrams': diagrams, 'checks': checks,
              'historical_dependency_hash_differences': [
                  {'orbit_id': d['orbit_id'], **dependency}
                  for d in diagrams for dependency in d['dependencies']
                  if dependency.get('declared_hash_matches') is False],
              'source_manifest': sorted(SOURCES.values(), key=lambda x: x['path'])}
    (output / 'diagram-ledger.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8', newline='\n')
    lines = ['# Three matter channels: per-diagram ledger', '',
             'This ledger contains 31 unique diagrams and 35 channel occurrences: 31 Af, two Ff and two ff. All coefficient positions from eps^-6 through eps^0 are explicit; source-metadata zeros remain distinguished from explicitly stored zeros.', '',
             '`Af = CF*CA*TF*nf`, `Ff = CF^2*TF*nf`, `ff = CF*(TF*nf)^2`.', '',
             'The expressions below are unit-scale physical web contributions. Multiply the channel sum by S3 once, then apply coupling counterterms. Native I already includes the adopted sign, C4 multiplicity, crossing and one Euler compensation. Pi1 and leg-Pi1 contributions require Bf(eps)=-8*(1-eps)/(3-2*eps); chain and Pi2 channel fields already contain their color operators.', '',
             'Here rho=log(z), r=1/(1+z), z>0. Raw diagram L means log(z); the renormalization-scale L belongs only to the subsequent S3/CT stage. H_z_WORD and H_r_WORD preserve the source Hlog argument and word, with m1 meaning -1. Their kernel is dt/(t-a); Hlog(r,[1])=log(1-r). Raw syntax, aliases, exact source definitions and compact alternatives are retained in [diagram-ledger.json](diagram-ledger.json).', '',
             'The 217 Af diagram-coefficient comparisons pass. The full mixed r-to-z Af sum uses the existing historical HyperInt witness; this export does not claim a fresh canonicalization, integration or formal proof.', '',
             '| Diagram | Branch | Channels | External multiplier |', '|---|---|---|---|']
    for d in diagrams:
        c = d['contributions']['Af']
        lines.append(f"| [{d['orbit_id']}](#{d['orbit_id'].lower().replace('+','')}) | {d['family']} | {', '.join(d['contributions'])} | `{c['external_scalar_multiplier']}` |")
    for d in diagrams:
        lines.extend(['', '## ' + d['orbit_id'], '',
                      f"Source: [{d['result_source']['path']}](sources/{d['result_source']['path']}) (revision {d['revision']}); SHA-256 `{d['result_source']['sha256']}`.", '',
                      f"Native unweighted field: `{d['native_unweighted']['source']['json_pointer']}`. {d['native_unweighted']['meaning']}.", '',
                      '### Native unweighted coefficients', ''])
        for q, value in d['native_unweighted']['coefficients'].items():
            lines.extend([f"`eps^{q}` — {value['origin']}:", '', '```text', value['normalized_expression'], '```', ''])
        for channel, c in d['contributions'].items():
            lines.extend(['### ' + channel + ' contribution', '',
                          f"Consumed field: `{c['consumed_native_field']['source']['json_pointer']}`. External scalar multiplier: `{c['external_scalar_multiplier']}`; physical assembly sign: +1.", ''])
            if c['color_operator_already_included']:
                lines.extend(['The consumed field already includes the channel color operator; do not multiply it again.', ''])
            for q, value in c['coefficients'].items():
                lines.extend([f'`eps^{q}`:', '', '```text', value['expression'], '```', ''])
        if d['compact_source_alternatives']:
            lines.extend(['Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.', ''])
        definitions = {k: v for k, v in d['source_metadata'].items()
                       if k in ['variables', 'Hlog_variables', 'Hlog_convention', 'functions',
                                'coefficient_notation', 'expression_language', 'Hlog_language',
                                'Hlog_expression_language', 'coefficient_field'] or k.endswith('_definition')}
        if definitions:
            lines.extend(['### Source-scoped symbol definitions', '',
                          'These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.', '',
                          '```json', json.dumps(definitions, ensure_ascii=False, indent=2), '```', ''])
        for dep in d['dependencies']:
            lines.append(f"Source dependency: [{dep['path']}](sources/{dep['path']}).")
        if any(dep.get('declared_hash_matches') is False for dep in d['dependencies']):
            lines.extend(['', 'A historical dependency digest differs from the current documentation file. Both digests are recorded in the JSON; no identity is inferred between the two versions. All selected coefficient-file digests match the fresh replay.', ''])
    (output / 'diagrams.md').write_text('\n'.join(lines).rstrip() + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({**result['summary'], 'source_files': len(SOURCES),
                      'output_dir': str(output)}, indent=2))


if __name__ == '__main__':
    main()
