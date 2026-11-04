# Complete adopted Pi1 moment evaluations

All functions below use fixed positive real z, real rho=log z, the HLOG kernel dt/(t-a), base point0 and endpointz. Write Hlog(z,[a1,...,an])=G(a1,...,an;z), zeta[n]=zeta(n), and use the empty word for constants. The evaluated functions are adopted existing mathematical results, not fresh integrations performed by the current Python replay.

## Ten original no-factorial moments

The exact adopted evaluation lemma is
\[
M_{k,j}(z)=\int_0^1\frac{\log[\beta(1-\beta)]^k[\rho+\log((1-\beta)/\beta)]^j}{\beta-z/(1+z)}d\beta
\]
with the ordinary removable value at the internal point and the ordinary integrable endpoint limits, for the ten index pairs printed below. Its complete values include all constants. The finite derivative rule(5) in [the main derivation](../pi1.md) follows by an actual linear endpoint primitive; substitution of these values evaluates every needed H jet through order3. In the original calculation, logarithmic products were expanded into hyperlogarithms, partial fractions supplied dt/(t-a) primitives, and the complete endpoint limits retained their zeta periods before fibration to z. The source [basis derivation](../runtime/snapshots/benchmarks/2ne/families/single-opposite/basis.md) and unchanged native files carry that adopted function evaluation. No unproved finite-part splitting at the removable internal point is introduced here.

The following are full native expressions; their notation has the exact definition above, and there is no missing word table:

### M_0,1

```text
M_0_1 = -3*zeta[2]-Hlog(z,[0, 0])
```

### M_0,2

```text
M_0_2 = -2*Hlog(z,[0, 0, 0])-2*zeta[2]*Hlog(z,[0])
```

### M_0,3

```text
M_0_3 = -9*zeta[2]^2-6*zeta[2]*Hlog(z,[0, 0])-6*Hlog(z,[0, 0, 0, 0])
```

### M_0,4

```text
M_0_4 = -24*zeta[2]*Hlog(z,[0, 0, 0])-24*Hlog(z,[0, 0, 0, 0, 0])-84/5*zeta[2]^2*Hlog(z,[0])
```

### M_1,1

```text
M_1_1 = -Hlog(z,[0, 0, 0])+2*zeta[3]+6*zeta[2]*Hlog(z,[-1])-3*zeta[2]*Hlog(z,[0])+2*Hlog(z,[-1, 0, 0])
```

### M_1,2

```text
M_1_2 = 4*zeta[2]*Hlog(z,[-1, 0])+17/5*zeta[2]^2-2*zeta[2]*Hlog(z,[0, 0])+4*zeta[3]*Hlog(z,[0])-2*Hlog(z,[0, 0, 0, 0])+4*Hlog(z,[-1, 0, 0, 0])
```

### M_1,3

```text
M_1_3 = -6*zeta[2]*Hlog(z,[0, 0, 0])+12*zeta[2]*zeta[3]+12*zeta[5]-6*Hlog(z,[0, 0, 0, 0, 0])+12*zeta[3]*Hlog(z,[0, 0])+12*Hlog(z,[-1, 0, 0, 0, 0])+18*zeta[2]^2*Hlog(z,[-1])-9*zeta[2]^2*Hlog(z,[0])+12*zeta[2]*Hlog(z,[-1, 0, 0])
```

### M_2,1

```text
M_2_1 = 12*zeta[2]*Hlog(z,[-1, 0])+4*Hlog(z,[0, -1, 0, 0])-24*zeta[2]*Hlog(z,[-1, -1])-8*Hlog(z,[-1, -1, 0, 0])-33/5*zeta[2]^2-4*zeta[2]*Hlog(z,[0, 0])-8*zeta[3]*Hlog(z,[-1])+4*zeta[3]*Hlog(z,[0])-2*Hlog(z,[0, 0, 0, 0])+4*Hlog(z,[-1, 0, 0, 0])+12*zeta[2]*Hlog(z,[0, -1])
```

### M_2,2

```text
M_2_2 = -16*zeta[3]*Hlog(z,[-1, 0])-16*zeta[2]*Hlog(z,[-1, -1, 0])+8*zeta[2]*zeta[3]+8*zeta[5]-16*Hlog(z,[-1, -1, 0, 0, 0])-4*Hlog(z,[0, 0, 0, 0, 0])+8*zeta[3]*Hlog(z,[0, 0])+8*Hlog(z,[0, -1, 0, 0, 0])+8*Hlog(z,[-1, 0, 0, 0, 0])-68/5*zeta[2]^2*Hlog(z,[-1])+6*zeta[2]^2*Hlog(z,[0])+8*zeta[2]*Hlog(z,[-1, 0, 0])+8*zeta[2]*Hlog(z,[0, -1, 0])
```

### M_3,1

```text
M_3_1 = -24*zeta[3]*Hlog(z,[-1, 0])+48*zeta[3]*Hlog(z,[-1, -1])-12*zeta[2]*Hlog(z,[0, 0, 0])-72*zeta[2]*Hlog(z,[0, -1, -1])+12*Hlog(z,[0, 0, -1, 0, 0])-72*zeta[2]*Hlog(z,[-1, -1, 0])-12*zeta[2]*zeta[3]+36*zeta[5]+48*Hlog(z,[-1, -1, -1, 0, 0])-24*Hlog(z,[0, -1, -1, 0, 0])+36*zeta[2]*Hlog(z,[0, 0, -1])+144*zeta[2]*Hlog(z,[-1, -1, -1])-24*Hlog(z,[-1, -1, 0, 0, 0])-6*Hlog(z,[0, 0, 0, 0, 0])+12*Hlog(z,[0, -1, 0, 0, 0])-72*zeta[2]*Hlog(z,[-1, 0, -1])+12*Hlog(z,[-1, 0, 0, 0, 0])+198/5*zeta[2]^2*Hlog(z,[-1])-24*Hlog(z,[-1, 0, -1, 0, 0])-99/5*zeta[2]^2*Hlog(z,[0])+24*zeta[2]*Hlog(z,[-1, 0, 0])+36*zeta[2]*Hlog(z,[0, -1, 0])-24*zeta[3]*Hlog(z,[0, -1])
```

## The complete first ordered double-kernel coefficient

Define the actual coupled body, with its own ordering retained,
\[
T(a,b;z)=(1+z)^2z^{-(a+b)/2}\int_{0<s_1<s_2<1\atop0<t_1<t_2<1}
D_z(s_1,t_1)^{a-1}D_z(s_2,t_2)^{b-1}d^4s,
\quad \rho_q=[e^q]T(e,e;z).
\]
Each opposite kernel has distinct endpoints; bounding by the full four-cube gives local holomorphy for small a,b. Its exact needed density is
\[
\rho_1=(1+z)^2\int_{\Omega_2}\frac{\log D_1+\log D_2-\log z}{D_1D_2}d^4s.
\]
The original outer-to-inner order is s2,s1,t2,t1; the corner bound justifies regrouping for the following actual reduction. Integrate t1 from0 to u=t2 first, put v=s1,w=s2 and set D1=D_z(u,v), D2=D_z(u,w), B=1-v, d=(1+z)v-1. The exact path primitive is (D1^e-B^e)/(e*d), with removable d=0 limit. Its Taylor expansion gives
\[
K_0=(1+z)^2\frac{\log(D_1/B)}{dD_2},\qquad
\rho_1=\int_0^1du\int_0^1dw\int_0^wdv\,K_0
[\log D_2+\tfrac12(\log D_1+\log B)-\log z].
\]
The adopted computational lemma equates this complete convergent kernel to the HLOG expression below, with ordinary combined endpoints and positive-real branches. The recorded finite integration route uses v,u,w after a=1/(1+z); the complete result is subsequently fibred to z. No separately singular term at d=0 is admitted. The original two-opposite HyperInt evaluation retains all endpoint constants. The [original method](../runtime/snapshots/benchmarks/2ne/families/same-channel-two/method.md) and [structural dictionary](../runtime/snapshots/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json) retain the exact evaluation source.

The solved zero coefficient is rho0=M01^2/6. For rho1 the full answer is

```text
rho1 = -12*zeta(5) + 17*pi**2*zeta(3)/6 - pi**4*Hlog(z, (-1,))/12 + 49*pi**4*Hlog(z, (0,))/360 + 3*zeta(3)*Hlog(z, (0, 0)) - pi**2*Hlog(z, (-1, 0, 0))/3 + pi**2*Hlog(z, (0, -1, 0))/6 - pi**2*Hlog(z, (0, 0, -1))/2 + pi**2*Hlog(z, (0, 0, 0))/3 - 2*Hlog(z, (-1, 0, 0, 0, 0)) + Hlog(z, (0, -1, 0, 0, 0)) - Hlog(z, (0, 0, -1, 0, 0)) + Hlog(z, (0, 0, 0, 0, 0))
```

The symbols rho0/rho1 here label Taylor coefficients and are distinct from the unindexed shape logarithm rho. At first order only,
\[
T(2e,e)+T(e,2e)=2\rho_0+3e\rho_1+O(e^2).
\]
The complete thirteen-term expression, not its dictionary name, is substituted into the actual P5 finite coefficient. This is an explicit adopted evaluated lemma; file hashes or an old successful CAS process alone would not establish the original integral.
