"""Display the existing compressed double-TGV source through its full five-line decoder."""
def emit_integral(orbit,info):
    from .definition_exporter import _build_orbit_definition
    from .emitter import _emit_measure,_emit_eikonal_chain
    raw=_build_orbit_definition(orbit)
    loops=raw['internal_vertices'][0]['loops']
    if len(loops)!=4:raise ValueError('double_vertex_requires_four_ordered_external_loops')
    by={p['loop']:p for p in raw['propagators']}
    phases={p['loop']:p['expression'] for p in raw['phases']}
    chain,_=_emit_eikonal_chain(info.edge_counts)
    text=[rf'I_{{({orbit.orbit_id})}} = '+_emit_measure(info.edge_counts)+
          r'\;\frac1{N_c}\operatorname{Tr}\bigl['+chain+r'\bigr]',
          r'\times\int\frac{d^Dq}{(2\pi)^D}\prod_{j=1}^4\int\frac{d^D\ell_j}{(2\pi)^D}'
          r'\;(2\pi)^{2D}\delta^{(D)}(\ell_1+\ell_2+q)\delta^{(D)}(\ell_3+\ell_4-q)',
          r'\times\mathcal V_3^{b_1b_2c;\nu_1\nu_2\rho}(-\ell_1,-\ell_2,-q)'
          r'\mathcal V_3^{b_3b_4d;\nu_3\nu_4\sigma}(-\ell_3,-\ell_4,q)'
          r'\frac{-ig_{\rho\sigma}\delta^{cd}}{q^2+i0}']
    for j,loop in enumerate(loops,1):
        prop=by[loop];mu=prop['lorentz'][0];a=prop['color'][0]
        def subscript(label):
            import re
            m=re.fullmatch(r'(mu|a)(\d+)',label)
            return ('\\mu' if m[1]=='mu' else 'a')+'_{'+m[2]+'}' if m else label
        text.append(r'\times\frac{-ig_{'+subscript(mu)+rf'\nu_{j}'+r'}\delta^{'+subscript(a)+rf'b_{j}'+'}}'
                    +rf'{{\ell_{j}^2+i0}}e^{{-i\ell_{j}\cdot\bigl('+phases[loop]+r'\bigr)}')
    return '\n'.join(text)
