import sympy as sp
from sympy import symbols, expand, simplify, collect, factor

# ============================================================
# Full Pnum computation for [1,1,1,1] (12)(34) via framework recipe
# ============================================================

z = symbols('z')
s1, s2, s3, s4 = symbols('s1 s2 s3 s4')
a1, a2, a3, a4, a5 = symbols('alpha1 alpha2 alpha3 alpha4 alpha5')
eps = symbols('eps')
D_val = 4 - 2*eps

# Kinematics
vd = {}
vd[(1,1)] = vd[(2,2)] = vd[(3,3)] = vd[(4,4)] = 0
vd[(1,2)] = vd[(2,1)] = z
vd[(3,4)] = vd[(4,3)] = z
vd[(1,4)] = vd[(4,1)] = 1
vd[(2,3)] = vd[(3,2)] = 1
vd[(1,3)] = vd[(3,1)] = -(1+z)
vd[(2,4)] = vd[(4,2)] = -(1+z)

def vdot(p, q):
    return vd[(p, q)]

# c_ij matrix
def cij(i, j):
    if i > j: i, j = j, i
    if (i,j) == (1,2): return a3 + a4 + a5
    elif (i,j) == (3,4): return a1 + a2 + a5
    elif i in (1,2) and j in (3,4): return a5
    return 0

# Edge of each leg
edge_of = {1:1, 2:2, 3:3, 4:4}

# V_i coefficients
V = {1: {1:0, 2:0, 3:0, 4:0},
     2: {1:1, 2:0, 3:0, 4:0},
     3: {1:1, 2:1, 3:0, 4:0},
     4: {1:1, 2:1, 3:1, 4:0}}

s_var = {1: s1, 2: s2, 3: s3, 4: s4}

def vdotX(p_edge, i_leg):
    e_i = edge_of[i_leg]
    si = s_var[e_i]
    result = 0
    for k in range(1, 5):
        result += V[e_i][k] * vdot(p_edge, k)
    result += si * vdot(p_edge, e_i)
    return expand(result)

def vdotXdiff(p_edge, i, j):
    return expand(vdotX(p_edge, i) - vdotX(p_edge, j))

alpha = {1: a1, 2: a2, 3: a3, 4: a4, 5: a5}

def F_ip(p_edge, i):
    ai = alpha[i]
    result = 0
    for j in range(1, 5):
        if j == i: continue
        aj = alpha[j]
        result += aj * cij(i, j) * vdotXdiff(p_edge, i, j)
    return -2 * ai * expand(result)

def H_ijpq(p_edge, q_edge, i, j):
    pq = vdot(p_edge, q_edge)
    if i != j:
        return 2 * alpha[i] * alpha[j] * cij(i, j) * pq
    else:
        s = 0
        for k in range(1, 5):
            if k == i: continue
            s += alpha[k] * cij(i, k)
        return -2 * alpha[i] * pq * s

def FG_typeI(i, j):
    """Compute (∂_i B_5)·(∂_j B_5)"""
    result = 0
    ai = alpha[i]
    aj = alpha[j]
    for k in range(1, 5):
        if k == i: continue
        ak = alpha[k]
        for l in range(1, 5):
            if l == j: continue
            al = alpha[l]
            # Build coefficient vector for X_i - X_k
            coeffs = {}
            for m in range(1, 5):
                coeff = V[edge_of[i]].get(m, 0) - V[edge_of[k]].get(m, 0)
                if m == edge_of[i]: coeff += s_var[edge_of[i]]
                if m == edge_of[k]: coeff -= s_var[edge_of[k]]
                if coeff != 0: coeffs[m] = coeff
            # Build coefficient vector for X_j - X_l
            coeffs2 = {}
            for n in range(1, 5):
                coeff = V[edge_of[j]].get(n, 0) - V[edge_of[l]].get(n, 0)
                if n == edge_of[j]: coeff += s_var[edge_of[j]]
                if n == edge_of[l]: coeff -= s_var[edge_of[l]]
                if coeff != 0: coeffs2[n] = coeff
            # Compute dot product
            dotprod = 0
            for m, cm in coeffs.items():
                for n, cn in coeffs2.items():
                    dotprod += cm * cn * vdot(m, n)
            result += ak * al * cij(i, k) * cij(j, l) * dotprod
    return 4 * ai * aj * expand(result)

def H_typeI(i, j):
    if i != j:
        return 2 * alpha[i] * alpha[j] * cij(i, j) * D_val
    else:
        s = 0
        for k in range(1, 5):
            if k == i: continue
            s += alpha[k] * cij(i, k)
        return -2 * alpha[i] * D_val * s

# ============================================================
# Operator term table
# ============================================================

terms = []

# Type I (4 terms): ('I', C, i, j)
terms.append(('I', -vdot(1,2)*vdot(3,4), 2, 4))
terms.append(('I', +vdot(1,2)*vdot(3,4), 2, 3))
terms.append(('I', +vdot(1,2)*vdot(3,4), 1, 4))
terms.append(('I', -vdot(1,2)*vdot(3,4), 1, 3))

# Type IIa (4 terms): ('v', C, p_edge, d_p, q_edge, d_q)
terms.append(('v', -2*vdot(1,2), 3, {1:1}, 4, {2:1}))
terms.append(('v', +2*vdot(1,2), 3, {2:1}, 4, {1:1}))
terms.append(('v', -vdot(1,2), 3, {2:1, 1:-1}, 4, {3:1}))
terms.append(('v', -vdot(1,2), 4, {1:1, 2:-1}, 3, {4:1}))

# Type IIb (2 terms)
terms.append(('v', -4*vdot(1,2), 4, {1:1}, 3, {2:1}))
terms.append(('v', +4*vdot(1,2), 3, {1:1}, 4, {2:1}))

# Type IIc (2 terms)
terms.append(('v', -vdot(3,4), 1, {1:1, 2:2}, 2, {3:1, 4:-1}))
terms.append(('v', +vdot(3,4), 2, {1:2, 2:1}, 1, {3:1, 4:-1}))

# Type IIIa (4 terms)
terms.append(('v', -2*vdot(1,3), 2, {1:2, 2:1}, 4, {1:1, 2:1}))
terms.append(('v', +2*vdot(1,4), 2, {1:2, 2:1}, 3, {1:1, 2:1}))
terms.append(('v', +2*vdot(2,3), 1, {1:1, 2:2}, 4, {1:1, 2:1}))
terms.append(('v', -2*vdot(2,4), 1, {1:1, 2:2}, 3, {1:1, 2:1}))

# Type IIIb (4 terms)
terms.append(('v', -vdot(1,3), 2, {1:2, 2:1}, 4, {3:1, 1:-1, 2:-1}))
terms.append(('v', -vdot(1,4), 2, {1:2, 2:1}, 3, {1:1, 2:1, 4:-1}))
terms.append(('v', +vdot(2,3), 1, {1:1, 2:2}, 4, {3:1, 1:-1, 2:-1}))
terms.append(('v', +vdot(2,4), 1, {1:1, 2:2}, 3, {1:1, 2:1, 4:-1}))

# ============================================================
# Compute Pnum^(1) and Pnum^(2)
# ============================================================

Pnum1 = 0
Pnum2 = 0

for t_idx, term in enumerate(terms):
    if term[0] == 'I':
        _, C, i, j = term
        C = expand(C)
        if C == 0: continue
        FG = FG_typeI(i, j)
        H = H_typeI(i, j)
        Pnum1 += C * FG
        Pnum2 += C * H
    else:
        _, C, p_edge, d_p, q_edge, d_q = term
        C = expand(C)
        if C == 0: continue
        Fp = sum(di * F_ip(p_edge, i) for i, di in d_p.items())
        Gq = sum(dj * F_ip(q_edge, j) for j, dj in d_q.items())
        Hpq = 0
        for i, di in d_p.items():
            for j, dj in d_q.items():
                Hpq += di * dj * H_ijpq(p_edge, q_edge, i, j)
        Pnum1 += C * expand(Fp * Gq)
        Pnum2 += C * expand(Hpq)

Pnum1 = expand(Pnum1)
Pnum2 = expand(Pnum2)

print('=== Pnum^(1) (before (2-3eps) factor) ===')
print(f'Total terms: {len(str(Pnum1))} chars')
# Show by alpha5 power
print()
for k in range(5):
    coeff = Pnum1.coeff(a5, k)
    if coeff != 0:
        c_exp = expand(coeff)
        print(f'alpha5^{k}: {c_exp}')

print()
print('=== Pnum^(2) (before minus sign) ===')
print(f'Total terms: {len(str(Pnum2))} chars')
for k in range(5):
    coeff = Pnum2.coeff(a5, k)
    if coeff != 0:
        c_exp = expand(coeff)
        print(f'alpha5^{k}: {c_exp}')