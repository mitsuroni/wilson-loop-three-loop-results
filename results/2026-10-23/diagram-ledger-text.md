# Three matter channels: per-diagram ledger

This ledger contains 31 unique diagrams and 35 channel occurrences: 31 Af, two Ff and two ff. All coefficient positions from eps^-6 through eps^0 are explicit; source-metadata zeros remain distinguished from explicitly stored zeros.

`Af = CF*CA*TF*nf`, `Ff = CF^2*TF*nf`, `ff = CF*(TF*nf)^2`.

The expressions below are unit-scale physical web contributions. Multiply the channel sum by S3 once, then apply coupling counterterms. Native I already includes the adopted sign, C4 multiplicity, crossing and one Euler compensation. Pi1 and leg-Pi1 contributions require Bf(eps)=-8*(1-eps)/(3-2*eps); chain and Pi2 channel fields already contain their color operators.

Here rho=log(z), r=1/(1+z), z>0. Raw diagram L means log(z); the renormalization-scale L belongs only to the subsequent S3/CT stage. H_z_WORD and H_r_WORD preserve the source Hlog argument and word, with m1 meaning -1. Their kernel is dt/(t-a); Hlog(r,[1])=log(1-r). Raw syntax, aliases, exact source definitions and compact alternatives are retained in [diagram-ledger.json](diagram-ledger.json).

The 217 Af diagram-coefficient comparisons pass. The full mixed r-to-z Af sum uses the existing historical HyperInt witness; this export does not claim a fresh canonicalization, integration or formal proof.

| Diagram | Branch | Channels | External multiplier |
|---|---|---|---|
| [4e-Pi1-01](#4e-pi1-01) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-02](#4e-pi1-02) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-03](#4e-pi1-03) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-04](#4e-pi1-04) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-05](#4e-pi1-05) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-06](#4e-pi1-06) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-07](#4e-pi1-07) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-08](#4e-pi1-08) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-09](#4e-pi1-09) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-10](#4e-pi1-10) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [4e-Pi1-11](#4e-pi1-11) | Pi1 | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [2e-Pi1+Pi1-01](#2e-pi1pi1-01) | chain | Af, ff | `1` |
| [2e-Pi1+Pi1-02](#2e-pi1pi1-02) | chain | Af, ff | `1` |
| [2e-Pi2-01](#2e-pi2-01) | Pi2 | Af, Ff | `1` |
| [2e-Pi2-02](#2e-pi2-02) | Pi2 | Af, Ff | `1` |
| [3e-3g+Pi1-01](#3e-3gpi1-01) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-02](#3e-3gpi1-02) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-03](#3e-3gpi1-03) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-04](#3e-3gpi1-04) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-05](#3e-3gpi1-05) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-06](#3e-3gpi1-06) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-07](#3e-3gpi1-07) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-08](#3e-3gpi1-08) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-09](#3e-3gpi1-09) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-10](#3e-3gpi1-10) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-11](#3e-3gpi1-11) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-3g+Pi1-12](#3e-3gpi1-12) | legs | Af | `(8*eps - 8)/(3 - 2*eps)` |
| [3e-Gamma1-01](#3e-gamma1-01) | Gamma | Af | `1` |
| [3e-Gamma1-02](#3e-gamma1-02) | Gamma | Af | `1` |
| [3e-Gamma1-03](#3e-gamma1-03) | Gamma | Af | `1` |
| [3e-Gamma1-04](#3e-gamma1-04) | Gamma | Af | `1` |

## 4e-Pi1-01

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-01/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-01/r002.json) (revision r002); SHA-256 `22da1c140fa1da316d4859388cf1d6fe66bcfcfc5df5e58c1b88629d71dcb0f5`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
-1/576
```

`eps^-4` — explicit source coefficient:

```text
-1/576
```

`eps^-3` — explicit source coefficient:

```text
-1/192 - pi**2/2304
```

`eps^-2` — explicit source coefficient:

```text
-5/576 - pi**2/2304 - zeta(3)/576
```

`eps^-1` — explicit source coefficient:

```text
-11/576 - pi**2/768 - 19*pi**4/276480 - zeta(3)/576
```

`eps^0` — explicit source coefficient:

```text
-7/192 - 5*pi**2/2304 - 19*pi**4/276480 - zeta(3)/192 - pi**2*zeta(3)/2304 - zeta(5)/960
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
1/216
```

`eps^-4`:

```text
1/324
```

`eps^-3`:

```text
11/972 + pi**2/864
```

`eps^-2`:

```text
zeta(3)/216 + pi**2/1296 + 49/2916
```

`eps^-1`:

```text
zeta(3)/324 + 19*pi**4/103680 + 11*pi**2/3888 + 341/8748
```

`eps^0`:

```text
zeta(5)/360 + 19*pi**4/155520 + 11*zeta(3)/972 + pi**2*zeta(3)/864 + 49*pi**2/11664 + 1897/26244
```

Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-01/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-01/derivation.md).

## 4e-Pi1-02

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-02/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-02/r002.json) (revision r002); SHA-256 `f4ef60739c56f6950b523a6b6755f5c19b4020f34f34564ca6480cc8198451c5`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
-1/576
```

`eps^-4` — explicit source coefficient:

```text
-1/576
```

`eps^-3` — explicit source coefficient:

```text
-1/192 - pi**2/2304
```

`eps^-2` — explicit source coefficient:

```text
-5/576 - pi**2/2304 - zeta(3)/576
```

`eps^-1` — explicit source coefficient:

```text
-11/576 - pi**2/768 - 19*pi**4/276480 - zeta(3)/576
```

`eps^0` — explicit source coefficient:

```text
-7/192 - 5*pi**2/2304 - 19*pi**4/276480 - zeta(3)/192 - pi**2*zeta(3)/2304 - zeta(5)/960
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
1/216
```

`eps^-4`:

```text
1/324
```

`eps^-3`:

```text
11/972 + pi**2/864
```

`eps^-2`:

```text
zeta(3)/216 + pi**2/1296 + 49/2916
```

`eps^-1`:

```text
zeta(3)/324 + 19*pi**4/103680 + 11*pi**2/3888 + 341/8748
```

`eps^0`:

```text
zeta(5)/360 + 19*pi**4/155520 + 11*zeta(3)/972 + pi**2*zeta(3)/864 + 49*pi**2/11664 + 1897/26244
```

Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-02/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-02/derivation.md).

## 4e-Pi1-03

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-03/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-03/r002.json) (revision r002); SHA-256 `928e181315e388c55416e507eadfc496fc07132eab9c656d80dd5ebb6b624354`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
-1/384
```

`eps^-3` — explicit source coefficient:

```text
-1/384
```

`eps^-2` — explicit source coefficient:

```text
pi**2*H_z_m1/64 + H_z_m1_0_0/32 - rho**3/384 - pi**2*rho/128 - 5*pi**2/1536 - 1/128 + zeta(3)/32
```

`eps^-1` — explicit source coefficient:

```text
5*pi**2*H_z_0_m1/128 + 5*H_z_0_m1_0_0/64 - 5*H_z_m1*zeta(3)/32 + pi**2*H_z_m1/32 + 5*pi**2*H_z_m1_0/128 + H_z_m1_0_0/16 + 5*H_z_m1_0_0_0/64 - 5*pi**2*H_z_m1_m1/64 - 5*H_z_m1_m1_0_0/32 - rho**4/1536 - rho**3/192 - 5*pi**2*rho**2/768 - pi**2*rho/64 + 5*rho*zeta(3)/64 - 7*pi**4/1536 - 5*pi**2/1536 - 5/384 + 41*zeta(3)/384
```

`eps^0` — explicit source coefficient:

```text
pi**2*H_z_0_0_m1/32 + H_z_0_0_m1_0_0/16 - 19*H_z_0_m1*zeta(3)/64 + 5*pi**2*H_z_0_m1/64 + pi**2*H_z_0_m1_0/32 + 5*H_z_0_m1_0_0/32 - 7*H_z_0_m1_0_0_0/64 - 19*pi**2*H_z_0_m1_m1/128 - 19*H_z_0_m1_m1_0_0/64 - 5*H_z_m1*zeta(3)/16 + pi**2*H_z_m1/16 + 167*pi**4*H_z_m1/7680 - 19*H_z_m1_0*zeta(3)/64 + 5*pi**2*H_z_m1_0/64 + H_z_m1_0_0/8 + pi**2*H_z_m1_0_0/32 + 5*H_z_m1_0_0_0/32 - 17*H_z_m1_0_0_0_0/64 - 19*pi**2*H_z_m1_0_m1/128 - 19*H_z_m1_0_m1_0_0/64 - 5*pi**2*H_z_m1_m1/32 + 19*H_z_m1_m1*zeta(3)/32 - 19*pi**2*H_z_m1_m1_0/128 - 5*H_z_m1_m1_0_0/16 - 19*H_z_m1_m1_0_0_0/64 + 19*pi**2*H_z_m1_m1_m1/64 + 19*H_z_m1_m1_m1_0_0/32 + 13*rho**5/5120 - rho**4/768 - rho**3/96 + 7*pi**2*rho**3/1536 - pi**2*rho**2/96 + rho**2*zeta(3)/128 - 167*pi**4*rho/15360 - pi**2*rho/32 + 5*rho*zeta(3)/32 - 133*pi**4/12288 - 5*pi**2/512 - 11/384 + 65*zeta(3)/384 + 35*zeta(5)/64
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
1/144
```

`eps^-3`:

```text
1/216
```

`eps^-2`:

```text
-pi**2*H_z_m1/24 - H_z_m1_0_0/12 + rho**3/144 + pi**2*rho/48 - zeta(3)/12 + 11/648 + 5*pi**2/576
```

`eps^-1`:

```text
-5*pi**2*H_z_0_m1/48 - 5*H_z_0_m1_0_0/24 - 5*pi**2*H_z_m1/72 + 5*H_z_m1*zeta(3)/12 - 5*pi**2*H_z_m1_0/48 - 5*H_z_m1_0_0/36 - 5*H_z_m1_0_0_0/24 + 5*pi**2*H_z_m1_m1/24 + 5*H_z_m1_m1_0_0/12 + rho**4/576 + 5*rho**3/432 + 5*pi**2*rho**2/288 - 5*rho*zeta(3)/24 + 5*pi**2*rho/144 - 37*zeta(3)/144 + 49/1944 + 5*pi**2/864 + 7*pi**4/576
```

`eps^0`:

```text
-pi**2*H_z_0_0_m1/12 - H_z_0_0_m1_0_0/6 - 25*pi**2*H_z_0_m1/144 + 19*H_z_0_m1*zeta(3)/24 - pi**2*H_z_0_m1_0/12 - 25*H_z_0_m1_0_0/72 + 7*H_z_0_m1_0_0_0/24 + 19*pi**2*H_z_0_m1_m1/48 + 19*H_z_0_m1_m1_0_0/24 - 167*pi**4*H_z_m1/2880 - 7*pi**2*H_z_m1/54 + 25*H_z_m1*zeta(3)/36 - 25*pi**2*H_z_m1_0/144 + 19*H_z_m1_0*zeta(3)/24 - pi**2*H_z_m1_0_0/12 - 7*H_z_m1_0_0/27 - 25*H_z_m1_0_0_0/72 + 17*H_z_m1_0_0_0_0/24 + 19*pi**2*H_z_m1_0_m1/48 + 19*H_z_m1_0_m1_0_0/24 - 19*H_z_m1_m1*zeta(3)/12 + 25*pi**2*H_z_m1_m1/72 + 19*pi**2*H_z_m1_m1_0/48 + 25*H_z_m1_m1_0_0/36 + 19*H_z_m1_m1_0_0_0/24 - 19*pi**2*H_z_m1_m1_m1/24 - 19*H_z_m1_m1_m1_0_0/12 - 13*rho**5/1920 + 5*rho**4/1728 - 7*pi**2*rho**3/576 + 7*rho**3/324 - rho**2*zeta(3)/48 + 19*pi**2*rho**2/864 - 25*rho*zeta(3)/72 + 7*pi**2*rho/108 + 167*pi**4*rho/5760 - 35*zeta(5)/24 - 73*zeta(3)/216 + 341/5832 + 55*pi**2/2592 + 343*pi**4/13824
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "Mkn": "existing M_{k,n}; definitions and evaluated files in dependencies"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis.md](sources/benchmarks/2ne/families/single-opposite/basis.md).
Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_2.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_2.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_4.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_4.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_2.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_2.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_3.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_2.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_2.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_3_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_3_1.txt).
Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-03/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-03/derivation.md).

## 4e-Pi1-04

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-04/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-04/r002.json) (revision r002); SHA-256 `cb67dd178a5f02cb9c3dbd9b34a29f47a4054fb0be19cb85d09687abcce16cd9`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
0
```

`eps^-2` — explicit source coefficient:

```text
pi**2*H_z_m1/128 + H_z_m1_0_0/64 - rho**3/768 - pi**2*rho/256 + zeta(3)/64
```

`eps^-1` — explicit source coefficient:

```text
pi**2*H_z_0_m1/64 + H_z_0_m1_0_0/32 - H_z_m1*zeta(3)/16 + pi**2*H_z_m1/128 + pi**2*H_z_m1_0/64 + H_z_m1_0_0/64 + H_z_m1_0_0_0/32 - pi**2*H_z_m1_m1/32 - H_z_m1_m1_0_0/16 + rho**4/3072 - rho**3/768 - pi**2*rho**2/384 - pi**2*rho/256 + rho*zeta(3)/32 - 37*pi**4/15360 + zeta(3)/64
```

`eps^0` — explicit source coefficient:

```text
pi**2*H_z_0_0_m1/256 + H_z_0_0_m1_0_0/128 - 13*H_z_0_m1*zeta(3)/128 + pi**2*H_z_0_m1/64 + pi**2*H_z_0_m1_0/256 + H_z_0_m1_0_0/32 - 5*H_z_0_m1_0_0_0/64 - 13*pi**2*H_z_0_m1_m1/256 - 13*H_z_0_m1_m1_0_0/128 - H_z_m1*zeta(3)/16 + 3*pi**2*H_z_m1/128 + 149*pi**4*H_z_m1/15360 - 13*H_z_m1_0*zeta(3)/128 + pi**2*H_z_m1_0/64 + 3*H_z_m1_0_0/64 + 3*pi**2*H_z_m1_0_0/256 + H_z_m1_0_0_0/32 - 17*H_z_m1_0_0_0_0/128 - 13*pi**2*H_z_m1_0_m1/256 - 13*H_z_m1_0_m1_0_0/128 - pi**2*H_z_m1_m1/32 + 13*H_z_m1_m1*zeta(3)/64 - 13*pi**2*H_z_m1_m1_0/256 - H_z_m1_m1_0_0/16 - 13*H_z_m1_m1_0_0_0/128 + 13*pi**2*H_z_m1_m1_m1/128 + 13*H_z_m1_m1_m1_0_0/64 + 13*rho**5/10240 + rho**4/3072 - rho**3/256 + pi**2*rho**3/384 - pi**2*rho**2/384 + rho**2*zeta(3)/64 - 149*pi**4*rho/30720 - 3*pi**2*rho/256 + rho*zeta(3)/32 - 37*pi**4/15360 + 3*zeta(3)/64 + pi**2*zeta(3)/128 + 29*zeta(5)/128
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
0
```

`eps^-2`:

```text
-pi**2*H_z_m1/48 - H_z_m1_0_0/24 + rho**3/288 + pi**2*rho/96 - zeta(3)/24
```

`eps^-1`:

```text
-pi**2*H_z_0_m1/24 - H_z_0_m1_0_0/12 - pi**2*H_z_m1/72 + H_z_m1*zeta(3)/6 - pi**2*H_z_m1_0/24 - H_z_m1_0_0/36 - H_z_m1_0_0_0/12 + pi**2*H_z_m1_m1/12 + H_z_m1_m1_0_0/6 - rho**4/1152 + rho**3/432 + pi**2*rho**2/144 - rho*zeta(3)/12 + pi**2*rho/144 - zeta(3)/36 + 37*pi**4/5760
```

`eps^0`:

```text
-pi**2*H_z_0_0_m1/96 - H_z_0_0_m1_0_0/48 - pi**2*H_z_0_m1/36 + 13*H_z_0_m1*zeta(3)/48 - pi**2*H_z_0_m1_0/96 - H_z_0_m1_0_0/18 + 5*H_z_0_m1_0_0_0/24 + 13*pi**2*H_z_0_m1_m1/96 + 13*H_z_0_m1_m1_0_0/48 - 149*pi**4*H_z_m1/5760 - 11*pi**2*H_z_m1/216 + H_z_m1*zeta(3)/9 - pi**2*H_z_m1_0/36 + 13*H_z_m1_0*zeta(3)/48 - pi**2*H_z_m1_0_0/32 - 11*H_z_m1_0_0/108 - H_z_m1_0_0_0/18 + 17*H_z_m1_0_0_0_0/48 + 13*pi**2*H_z_m1_0_m1/96 + 13*H_z_m1_0_m1_0_0/48 - 13*H_z_m1_m1*zeta(3)/24 + pi**2*H_z_m1_m1/18 + 13*pi**2*H_z_m1_m1_0/96 + H_z_m1_m1_0_0/9 + 13*H_z_m1_m1_0_0_0/48 - 13*pi**2*H_z_m1_m1_m1/48 - 13*H_z_m1_m1_m1_0_0/24 - 13*rho**5/3840 - rho**4/1728 - pi**2*rho**3/144 + 11*rho**3/1296 - rho**2*zeta(3)/24 + pi**2*rho**2/216 - rho*zeta(3)/18 + 11*pi**2*rho/432 + 149*pi**4*rho/11520 - 29*zeta(5)/48 - pi**2*zeta(3)/48 - 11*zeta(3)/108 + 37*pi**4/8640
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "tau0qk": "[eps^k]H0",
    "tau2qk": "[eps^k](H2(z)+H2(1/z))/2",
    "eta2qk": "log(z)*[eps^k](H2(z)-H2(1/z))/2"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [benchmarks/2ne/families/single-opposite/basis.md](sources/benchmarks/2ne/families/single-opposite/basis.md).
Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-04/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-04/derivation.md).

## 4e-Pi1-05

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-05/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-05/r002.json) (revision r002); SHA-256 `6869b5e35cb56aaf77f19b80ac55bab5e45869f06da10c705caf26e83494a2b4`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
-pi**2/384
```

`eps^-2` — explicit source coefficient:

```text
-pi**2/384 + 3*zeta(3)/64
```

`eps^-1` — explicit source coefficient:

```text
pi**2*rho**2/384 - 37*pi**4/23040 - pi**2/128 + 3*zeta(3)/64
```

`eps^0` — explicit source coefficient:

```text
-3*rho**2*zeta(3)/64 + pi**2*rho**2/384 - 37*pi**4/23040 - 5*pi**2/384 - 5*pi**2*zeta(3)/768 + 9*zeta(3)/64 + 21*zeta(5)/64
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
pi**2/144
```

`eps^-2`:

```text
-zeta(3)/8 + pi**2/216
```

`eps^-1`:

```text
-pi**2*rho**2/144 - zeta(3)/12 + 11*pi**2/648 + 37*pi**4/8640
```

`eps^0`:

```text
-pi**2*rho**2/216 + rho**2*zeta(3)/8 - 7*zeta(5)/8 - 11*zeta(3)/36 + 5*pi**2*zeta(3)/288 + 49*pi**2/1944 + 37*pi**4/12960
```

Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-05/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-05/derivation.md).

## 4e-Pi1-06

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-06/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-06/r002.json) (revision r002); SHA-256 `9da31233ebc761d2d659ec28ab4e8022c42f730865d49e8d2cf3abe66fa084f6`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
-pi**2/384
```

`eps^-2` — explicit source coefficient:

```text
-pi**2/384 + 3*zeta(3)/64
```

`eps^-1` — explicit source coefficient:

```text
pi**2*rho**2/384 - 37*pi**4/23040 - pi**2/128 + 3*zeta(3)/64
```

`eps^0` — explicit source coefficient:

```text
-3*rho**2*zeta(3)/64 + pi**2*rho**2/384 - 37*pi**4/23040 - 5*pi**2/384 - 5*pi**2*zeta(3)/768 + 9*zeta(3)/64 + 21*zeta(5)/64
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
pi**2/144
```

`eps^-2`:

```text
-zeta(3)/8 + pi**2/216
```

`eps^-1`:

```text
-pi**2*rho**2/144 - zeta(3)/12 + 11*pi**2/648 + 37*pi**4/8640
```

`eps^0`:

```text
-pi**2*rho**2/216 + rho**2*zeta(3)/8 - 7*zeta(5)/8 - 11*zeta(3)/36 + 5*pi**2*zeta(3)/288 + 49*pi**2/1944 + 37*pi**4/12960
```

Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-06/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-06/derivation.md).

## 4e-Pi1-07

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-07/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-07/r002.json) (revision r002); SHA-256 `baafb42f20e8266ad5b5eca46db2d00eac1c2727226e8c5dd86eb1812517d8d9`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
0
```

`eps^-2` — explicit source coefficient:

```text
0
```

`eps^-1` — explicit source coefficient:

```text
pi**2*H_z_m1/256 + H_z_m1_0_0/128 - rho**4/1536 - rho**3/1536 - pi**2*rho**2/768 - pi**2*rho/512 - pi**4/1536 + zeta(3)/128
```

`eps^0` — explicit source coefficient:

```text
3*pi**2*H_z_0_0_m1/256 + 3*H_z_0_0_m1_0_0/128 + pi**2*H_z_0_m1/128 - pi**2*H_z_0_m1_0/256 + H_z_0_m1_0_0/64 - 3*H_z_0_m1_0_0_0/128 - H_z_m1*zeta(3)/32 + pi**2*H_z_m1/256 + pi**4*H_z_m1/512 + pi**2*H_z_m1_0/128 + H_z_m1_0_0/128 + pi**2*H_z_m1_0_0/128 + H_z_m1_0_0_0/64 + 3*H_z_m1_0_0_0_0/64 - pi**2*H_z_m1_m1/64 - H_z_m1_m1_0_0/32 - rho**5/5120 - rho**4/768 - pi**2*rho**3/768 - rho**3/1536 - 9*rho**2*zeta(3)/256 - 13*pi**2*rho**2/3072 - 49*pi**4*rho/15360 - pi**2*rho/512 + rho*zeta(3)/64 - 17*pi**2*zeta(3)/256 - 41*pi**4/15360 + zeta(3)/128 + 9*zeta(5)/32
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
0
```

`eps^-2`:

```text
0
```

`eps^-1`:

```text
-pi**2*H_z_m1/96 - H_z_m1_0_0/48 + rho**4/576 + rho**3/576 + pi**2*rho**2/288 + pi**2*rho/192 - zeta(3)/48 + pi**4/576
```

`eps^0`:

```text
-pi**2*H_z_0_0_m1/32 - H_z_0_0_m1_0_0/16 - pi**2*H_z_0_m1/48 + pi**2*H_z_0_m1_0/96 - H_z_0_m1_0_0/24 + H_z_0_m1_0_0_0/16 - pi**4*H_z_m1/192 - pi**2*H_z_m1/144 + H_z_m1*zeta(3)/12 - pi**2*H_z_m1_0/48 - pi**2*H_z_m1_0_0/48 - H_z_m1_0_0/72 - H_z_m1_0_0_0/24 - H_z_m1_0_0_0_0/8 + pi**2*H_z_m1_m1/24 + H_z_m1_m1_0_0/12 + rho**5/1920 + 5*rho**4/1728 + rho**3/864 + pi**2*rho**3/288 + 35*pi**2*rho**2/3456 + 3*rho**2*zeta(3)/32 - rho*zeta(3)/24 + pi**2*rho/288 + 49*pi**4*rho/5760 - 3*zeta(5)/4 - zeta(3)/72 + 113*pi**4/17280 + 17*pi**2*zeta(3)/96
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "rho1": "existing rho_id_0_0_1 = [eps] Q_id;0,0(eps,z), explicit Hlog in dependency",
    "L": "log(z)"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [benchmarks/2ne/families/same-channel-two/structural-basis-functions.json](sources/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json).
Source dependency: [benchmarks/2ne/families/same-channel-two/unified-structure.md](sources/benchmarks/2ne/families/same-channel-two/unified-structure.md).
Source dependency: [benchmarks/2ne/families/same-channel-two/method.md](sources/benchmarks/2ne/families/same-channel-two/method.md).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-07/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-07/derivation.md).

## 4e-Pi1-08

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-08/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-08/r002.json) (revision r002); SHA-256 `c14243865814c45343f68cc8990f8793356c7c23d88e7a62058b69b599c100dd`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
0
```

`eps^-2` — explicit source coefficient:

```text
0
```

`eps^-1` — explicit source coefficient:

```text
pi**2*H_z_m1/256 + H_z_m1_0_0/128 - rho**4/1536 - rho**3/1536 - pi**2*rho**2/768 - pi**2*rho/512 - pi**4/1536 + zeta(3)/128
```

`eps^0` — explicit source coefficient:

```text
3*pi**2*H_z_0_0_m1/256 + 3*H_z_0_0_m1_0_0/128 + pi**2*H_z_0_m1/128 - pi**2*H_z_0_m1_0/256 + H_z_0_m1_0_0/64 - 3*H_z_0_m1_0_0_0/128 - H_z_m1*zeta(3)/32 + pi**2*H_z_m1/256 + pi**4*H_z_m1/512 + pi**2*H_z_m1_0/128 + H_z_m1_0_0/128 + pi**2*H_z_m1_0_0/128 + H_z_m1_0_0_0/64 + 3*H_z_m1_0_0_0_0/64 - pi**2*H_z_m1_m1/64 - H_z_m1_m1_0_0/32 - rho**5/5120 - rho**4/768 - pi**2*rho**3/768 - rho**3/1536 - 9*rho**2*zeta(3)/256 - 13*pi**2*rho**2/3072 - 49*pi**4*rho/15360 - pi**2*rho/512 + rho*zeta(3)/64 - 17*pi**2*zeta(3)/256 - 41*pi**4/15360 + zeta(3)/128 + 9*zeta(5)/32
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
0
```

`eps^-2`:

```text
0
```

`eps^-1`:

```text
-pi**2*H_z_m1/96 - H_z_m1_0_0/48 + rho**4/576 + rho**3/576 + pi**2*rho**2/288 + pi**2*rho/192 - zeta(3)/48 + pi**4/576
```

`eps^0`:

```text
-pi**2*H_z_0_0_m1/32 - H_z_0_0_m1_0_0/16 - pi**2*H_z_0_m1/48 + pi**2*H_z_0_m1_0/96 - H_z_0_m1_0_0/24 + H_z_0_m1_0_0_0/16 - pi**4*H_z_m1/192 - pi**2*H_z_m1/144 + H_z_m1*zeta(3)/12 - pi**2*H_z_m1_0/48 - pi**2*H_z_m1_0_0/48 - H_z_m1_0_0/72 - H_z_m1_0_0_0/24 - H_z_m1_0_0_0_0/8 + pi**2*H_z_m1_m1/24 + H_z_m1_m1_0_0/12 + rho**5/1920 + 5*rho**4/1728 + rho**3/864 + pi**2*rho**3/288 + 35*pi**2*rho**2/3456 + 3*rho**2*zeta(3)/32 - rho*zeta(3)/24 + pi**2*rho/288 + 49*pi**4*rho/5760 - 3*zeta(5)/4 - zeta(3)/72 + 113*pi**4/17280 + 17*pi**2*zeta(3)/96
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "rho1": "existing rho_id_0_0_1 = [eps] Q_id;0,0(eps,z), explicit Hlog in dependency",
    "L": "log(z)"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [benchmarks/2ne/families/same-channel-two/structural-basis-functions.json](sources/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json).
Source dependency: [benchmarks/2ne/families/same-channel-two/unified-structure.md](sources/benchmarks/2ne/families/same-channel-two/unified-structure.md).
Source dependency: [benchmarks/2ne/families/same-channel-two/method.md](sources/benchmarks/2ne/families/same-channel-two/method.md).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-08/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-08/derivation.md).

## 4e-Pi1-09

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-09/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-09/r002.json) (revision r002); SHA-256 `753efcab2092a5954ce44bc0b91f8880fc9e1943b7fed9648b75085165ad85d3`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
0
```

`eps^-2` — explicit source coefficient:

```text
pi**2*H_z_m1/128 + H_z_m1_0_0/64 - rho**3/768 - pi**2*rho/256 + zeta(3)/64
```

`eps^-1` — explicit source coefficient:

```text
pi**2*H_z_0_m1/64 + H_z_0_m1_0_0/32 - H_z_m1*zeta(3)/16 + pi**2*H_z_m1/128 + pi**2*H_z_m1_0/64 + H_z_m1_0_0/64 + H_z_m1_0_0_0/32 - pi**2*H_z_m1_m1/32 - H_z_m1_m1_0_0/16 + rho**4/3072 - rho**3/768 - pi**2*rho**2/384 - pi**2*rho/256 + rho*zeta(3)/32 - 37*pi**4/15360 + zeta(3)/64
```

`eps^0` — explicit source coefficient:

```text
pi**2*H_z_0_0_m1/256 + H_z_0_0_m1_0_0/128 - 13*H_z_0_m1*zeta(3)/128 + pi**2*H_z_0_m1/64 + pi**2*H_z_0_m1_0/256 + H_z_0_m1_0_0/32 - 5*H_z_0_m1_0_0_0/64 - 13*pi**2*H_z_0_m1_m1/256 - 13*H_z_0_m1_m1_0_0/128 - H_z_m1*zeta(3)/16 + 3*pi**2*H_z_m1/128 + 149*pi**4*H_z_m1/15360 - 13*H_z_m1_0*zeta(3)/128 + pi**2*H_z_m1_0/64 + 3*H_z_m1_0_0/64 + 3*pi**2*H_z_m1_0_0/256 + H_z_m1_0_0_0/32 - 17*H_z_m1_0_0_0_0/128 - 13*pi**2*H_z_m1_0_m1/256 - 13*H_z_m1_0_m1_0_0/128 - pi**2*H_z_m1_m1/32 + 13*H_z_m1_m1*zeta(3)/64 - 13*pi**2*H_z_m1_m1_0/256 - H_z_m1_m1_0_0/16 - 13*H_z_m1_m1_0_0_0/128 + 13*pi**2*H_z_m1_m1_m1/128 + 13*H_z_m1_m1_m1_0_0/64 + 13*rho**5/10240 + rho**4/3072 - rho**3/256 + pi**2*rho**3/384 - pi**2*rho**2/384 + rho**2*zeta(3)/64 - 149*pi**4*rho/30720 - 3*pi**2*rho/256 + rho*zeta(3)/32 - 37*pi**4/15360 + 3*zeta(3)/64 + pi**2*zeta(3)/128 + 29*zeta(5)/128
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
0
```

`eps^-2`:

```text
-pi**2*H_z_m1/48 - H_z_m1_0_0/24 + rho**3/288 + pi**2*rho/96 - zeta(3)/24
```

`eps^-1`:

```text
-pi**2*H_z_0_m1/24 - H_z_0_m1_0_0/12 - pi**2*H_z_m1/72 + H_z_m1*zeta(3)/6 - pi**2*H_z_m1_0/24 - H_z_m1_0_0/36 - H_z_m1_0_0_0/12 + pi**2*H_z_m1_m1/12 + H_z_m1_m1_0_0/6 - rho**4/1152 + rho**3/432 + pi**2*rho**2/144 - rho*zeta(3)/12 + pi**2*rho/144 - zeta(3)/36 + 37*pi**4/5760
```

`eps^0`:

```text
-pi**2*H_z_0_0_m1/96 - H_z_0_0_m1_0_0/48 - pi**2*H_z_0_m1/36 + 13*H_z_0_m1*zeta(3)/48 - pi**2*H_z_0_m1_0/96 - H_z_0_m1_0_0/18 + 5*H_z_0_m1_0_0_0/24 + 13*pi**2*H_z_0_m1_m1/96 + 13*H_z_0_m1_m1_0_0/48 - 149*pi**4*H_z_m1/5760 - 11*pi**2*H_z_m1/216 + H_z_m1*zeta(3)/9 - pi**2*H_z_m1_0/36 + 13*H_z_m1_0*zeta(3)/48 - pi**2*H_z_m1_0_0/32 - 11*H_z_m1_0_0/108 - H_z_m1_0_0_0/18 + 17*H_z_m1_0_0_0_0/48 + 13*pi**2*H_z_m1_0_m1/96 + 13*H_z_m1_0_m1_0_0/48 - 13*H_z_m1_m1*zeta(3)/24 + pi**2*H_z_m1_m1/18 + 13*pi**2*H_z_m1_m1_0/96 + H_z_m1_m1_0_0/9 + 13*H_z_m1_m1_0_0_0/48 - 13*pi**2*H_z_m1_m1_m1/48 - 13*H_z_m1_m1_m1_0_0/24 - 13*rho**5/3840 - rho**4/1728 - pi**2*rho**3/144 + 11*rho**3/1296 - rho**2*zeta(3)/24 + pi**2*rho**2/216 - rho*zeta(3)/18 + 11*pi**2*rho/432 + 149*pi**4*rho/11520 - 29*zeta(5)/48 - pi**2*zeta(3)/48 - 11*zeta(3)/108 + 37*pi**4/8640
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "tau0qk": "[eps^k]H0",
    "tau2qk": "[eps^k](H2(z)+H2(1/z))/2",
    "eta2qk": "log(z)*[eps^k](H2(z)-H2(1/z))/2"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [benchmarks/2ne/families/single-opposite/basis.md](sources/benchmarks/2ne/families/single-opposite/basis.md).
Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-09/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-09/derivation.md).

## 4e-Pi1-10

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-10/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-10/r002.json) (revision r002); SHA-256 `75654d49701d92ad1c0fca96febd12d6ee9d8ae296a84525aa80ae89dfc2ebe4`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
-1/384
```

`eps^-3` — explicit source coefficient:

```text
-1/384
```

`eps^-2` — explicit source coefficient:

```text
pi**2*H_z_m1/64 + H_z_m1_0_0/32 - rho**3/384 - pi**2*rho/128 - 5*pi**2/1536 - 1/128 + zeta(3)/32
```

`eps^-1` — explicit source coefficient:

```text
5*pi**2*H_z_0_m1/128 + 5*H_z_0_m1_0_0/64 - 5*H_z_m1*zeta(3)/32 + pi**2*H_z_m1/32 + 5*pi**2*H_z_m1_0/128 + H_z_m1_0_0/16 + 5*H_z_m1_0_0_0/64 - 5*pi**2*H_z_m1_m1/64 - 5*H_z_m1_m1_0_0/32 - rho**4/1536 - rho**3/192 - 5*pi**2*rho**2/768 - pi**2*rho/64 + 5*rho*zeta(3)/64 - 7*pi**4/1536 - 5*pi**2/1536 - 5/384 + 41*zeta(3)/384
```

`eps^0` — explicit source coefficient:

```text
pi**2*H_z_0_0_m1/32 + H_z_0_0_m1_0_0/16 - 19*H_z_0_m1*zeta(3)/64 + 5*pi**2*H_z_0_m1/64 + pi**2*H_z_0_m1_0/32 + 5*H_z_0_m1_0_0/32 - 7*H_z_0_m1_0_0_0/64 - 19*pi**2*H_z_0_m1_m1/128 - 19*H_z_0_m1_m1_0_0/64 - 5*H_z_m1*zeta(3)/16 + pi**2*H_z_m1/16 + 167*pi**4*H_z_m1/7680 - 19*H_z_m1_0*zeta(3)/64 + 5*pi**2*H_z_m1_0/64 + H_z_m1_0_0/8 + pi**2*H_z_m1_0_0/32 + 5*H_z_m1_0_0_0/32 - 17*H_z_m1_0_0_0_0/64 - 19*pi**2*H_z_m1_0_m1/128 - 19*H_z_m1_0_m1_0_0/64 - 5*pi**2*H_z_m1_m1/32 + 19*H_z_m1_m1*zeta(3)/32 - 19*pi**2*H_z_m1_m1_0/128 - 5*H_z_m1_m1_0_0/16 - 19*H_z_m1_m1_0_0_0/64 + 19*pi**2*H_z_m1_m1_m1/64 + 19*H_z_m1_m1_m1_0_0/32 + 13*rho**5/5120 - rho**4/768 - rho**3/96 + 7*pi**2*rho**3/1536 - pi**2*rho**2/96 + rho**2*zeta(3)/128 - 167*pi**4*rho/15360 - pi**2*rho/32 + 5*rho*zeta(3)/32 - 133*pi**4/12288 - 5*pi**2/512 - 11/384 + 65*zeta(3)/384 + 35*zeta(5)/64
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
1/144
```

`eps^-3`:

```text
1/216
```

`eps^-2`:

```text
-pi**2*H_z_m1/24 - H_z_m1_0_0/12 + rho**3/144 + pi**2*rho/48 - zeta(3)/12 + 11/648 + 5*pi**2/576
```

`eps^-1`:

```text
-5*pi**2*H_z_0_m1/48 - 5*H_z_0_m1_0_0/24 - 5*pi**2*H_z_m1/72 + 5*H_z_m1*zeta(3)/12 - 5*pi**2*H_z_m1_0/48 - 5*H_z_m1_0_0/36 - 5*H_z_m1_0_0_0/24 + 5*pi**2*H_z_m1_m1/24 + 5*H_z_m1_m1_0_0/12 + rho**4/576 + 5*rho**3/432 + 5*pi**2*rho**2/288 - 5*rho*zeta(3)/24 + 5*pi**2*rho/144 - 37*zeta(3)/144 + 49/1944 + 5*pi**2/864 + 7*pi**4/576
```

`eps^0`:

```text
-pi**2*H_z_0_0_m1/12 - H_z_0_0_m1_0_0/6 - 25*pi**2*H_z_0_m1/144 + 19*H_z_0_m1*zeta(3)/24 - pi**2*H_z_0_m1_0/12 - 25*H_z_0_m1_0_0/72 + 7*H_z_0_m1_0_0_0/24 + 19*pi**2*H_z_0_m1_m1/48 + 19*H_z_0_m1_m1_0_0/24 - 167*pi**4*H_z_m1/2880 - 7*pi**2*H_z_m1/54 + 25*H_z_m1*zeta(3)/36 - 25*pi**2*H_z_m1_0/144 + 19*H_z_m1_0*zeta(3)/24 - pi**2*H_z_m1_0_0/12 - 7*H_z_m1_0_0/27 - 25*H_z_m1_0_0_0/72 + 17*H_z_m1_0_0_0_0/24 + 19*pi**2*H_z_m1_0_m1/48 + 19*H_z_m1_0_m1_0_0/24 - 19*H_z_m1_m1*zeta(3)/12 + 25*pi**2*H_z_m1_m1/72 + 19*pi**2*H_z_m1_m1_0/48 + 25*H_z_m1_m1_0_0/36 + 19*H_z_m1_m1_0_0_0/24 - 19*pi**2*H_z_m1_m1_m1/24 - 19*H_z_m1_m1_m1_0_0/12 - 13*rho**5/1920 + 5*rho**4/1728 - 7*pi**2*rho**3/576 + 7*rho**3/324 - rho**2*zeta(3)/48 + 19*pi**2*rho**2/864 - 25*rho*zeta(3)/72 + 7*pi**2*rho/108 + 167*pi**4*rho/5760 - 35*zeta(5)/24 - 73*zeta(3)/216 + 341/5832 + 55*pi**2/2592 + 343*pi**4/13824
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "Mkn": "existing M_{k,n}; definitions and evaluated files in dependencies"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis.md](sources/benchmarks/2ne/families/single-opposite/basis.md).
Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_2.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_2.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_4.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_4.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_2.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_2.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_3.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_2.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_2.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_3_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_3_1.txt).
Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-10/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-10/derivation.md).

## 4e-Pi1-11

Source: [benchmarks/self-energy/pi1/results/4e-Pi1-11/r002.json](sources/benchmarks/self-energy/pi1/results/4e-Pi1-11/r002.json) (revision r002); SHA-256 `f04657a1aab9b963470d66587a2886628d2cb9bf9c59fa92f64a22f4fa666fef`.

Native unweighted field: `/laurent/coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
0
```

`eps^-2` — explicit source coefficient:

```text
rho**2/256 + pi**2/256
```

`eps^-1` — explicit source coefficient:

```text
-pi**2*H_z_m1/256 - H_z_m1_0_0/128 - rho**4/256 + rho**3/1536 - pi**2*rho**2/128 + rho**2/256 + pi**2*rho/512 - pi**4/256 - zeta(3)/128 + pi**2/256
```

`eps^0` — explicit source coefficient:

```text
3*pi**2*H_z_0_0_m1/128 + 3*H_z_0_0_m1_0_0/64 - pi**2*H_z_0_m1/512 + 3*pi**2*H_z_0_m1_0/128 - H_z_0_m1_0_0/256 + 9*H_z_0_m1_0_0_0/64 - pi**2*H_z_m1/256 + H_z_m1*zeta(3)/128 + 3*pi**4*H_z_m1/256 - pi**2*H_z_m1_0/512 - H_z_m1_0_0/128 + 3*pi**2*H_z_m1_0_0/64 - H_z_m1_0_0_0/256 + 9*H_z_m1_0_0_0_0/32 + pi**2*H_z_m1_m1/256 + H_z_m1_m1_0_0/128 - rho**5/512 - 31*rho**4/3072 - pi**2*rho**3/128 + rho**3/1536 - 17*pi**2*rho**2/1024 + 3*rho**2/256 + 3*rho**2*zeta(3)/128 - 3*pi**4*rho/512 - rho*zeta(3)/256 + pi**2*rho/512 - 101*pi**4/15360 - zeta(3)/128 + 3*pi**2/256 + 3*pi**2*zeta(3)/128
```

### Af contribution

Consumed field: `/laurent/coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
0
```

`eps^-2`:

```text
-rho**2/96 - pi**2/96
```

`eps^-1`:

```text
pi**2*H_z_m1/96 + H_z_m1_0_0/48 + rho**4/96 - rho**3/576 - rho**2/144 + pi**2*rho**2/48 - pi**2*rho/192 - pi**2/144 + zeta(3)/48 + pi**4/96
```

`eps^0`:

```text
-pi**2*H_z_0_0_m1/16 - H_z_0_0_m1_0_0/8 + pi**2*H_z_0_m1/192 - pi**2*H_z_0_m1_0/16 + H_z_0_m1_0_0/96 - 3*H_z_0_m1_0_0_0/8 - pi**4*H_z_m1/32 - H_z_m1*zeta(3)/48 + pi**2*H_z_m1/144 + pi**2*H_z_m1_0/192 - pi**2*H_z_m1_0_0/8 + H_z_m1_0_0/72 + H_z_m1_0_0_0/96 - 3*H_z_m1_0_0_0_0/4 - pi**2*H_z_m1_m1/96 - H_z_m1_m1_0_0/48 + rho**5/192 + 3*rho**4/128 - rho**3/864 + pi**2*rho**3/48 - rho**2*zeta(3)/16 - 11*rho**2/432 + 43*pi**2*rho**2/1152 - pi**2*rho/288 + rho*zeta(3)/96 + pi**4*rho/64 - pi**2*zeta(3)/16 - 11*pi**2/432 + zeta(3)/72 + 9*pi**4/640
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "T0": "tau0q0 = M01",
    "T1": "tau0q1 = M11/2",
    "T2": "tau0q2 = [delta^2]H0(delta,z)",
    "L": "log(z)",
    "H0": "z^(-1/2*delta)*(-z - 1)*MeijerG[{{0, 0, delta}, {}}, {{0, delta, delta}, {}}, z]/(Gamma[1 - delta]*Gamma[delta + 1]^2)"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [benchmarks/2ne/families/single-opposite/basis.md](sources/benchmarks/2ne/families/single-opposite/basis.md).
Source dependency: [benchmarks/self-energy/pi1/results/2e-Pi1-02/r001.json](sources/benchmarks/self-energy/pi1/results/2e-Pi1-02/r001.json).
Source dependency: [benchmarks/self-energy/pi1/workspaces/4e-Pi1-11/derivation.md](sources/benchmarks/self-energy/pi1/workspaces/4e-Pi1-11/derivation.md).

## 2e-Pi1+Pi1-01

Source: [benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-01/r002.json](sources/benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-01/r002.json) (revision r002); SHA-256 `1809dd8f94191aadcdf1ae48824324be8121a9bb5a083913ba7f1065b863e926`.

Native unweighted field: `/laurent/coefficients`. unweighted scalar I; chain color multipliers are already included in the consumed channel fields.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
-1/576
```

`eps^-3` — explicit source coefficient:

```text
-5/1152
```

`eps^-2` — explicit source coefficient:

```text
-1/64 - pi**2/2304
```

`eps^-1` — explicit source coefficient:

```text
-5/144 - 7*zeta(3)/576 - 5*pi**2/4608
```

`eps^0` — explicit source coefficient:

```text
-7/72 - pi**2/256 - 35*zeta(3)/1152 - 67*pi**4/276480
```

### Af contribution

Consumed field: `/color_channels/1/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
5/162
```

`eps^-3`:

```text
67/972
```

`eps^-2`:

```text
5*pi**2/648 + 61/243
```

`eps^-1`:

```text
67*pi**2/3888 + 35*zeta(3)/162 + 4577/8748
```

`eps^0`:

```text
67*pi**4/15552 + 469*zeta(3)/972 + 61*pi**2/972 + 19597/13122
```

### ff contribution

Consumed field: `/color_channels/2/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
-1/81
```

`eps^-3`:

```text
-11/486
```

`eps^-2`:

```text
-7/81 - pi**2/324
```

`eps^-1`:

```text
-703/4374 - 7*zeta(3)/81 - 11*pi**2/1944
```

`eps^0`:

```text
-3179/6561 - 7*pi**2/324 - 77*zeta(3)/486 - 67*pi**4/38880
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "Mkn": "existing evaluated M_{k,n}; bound Hlog files in dependencies",
    "H0": "z^(-1/2*delta)*(-z - 1)*MeijerG[{{0, 0, delta}, {}}, {{0, delta, delta}, {}}, z]/(Gamma[1 - delta]*Gamma[delta + 1]^2)"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [benchmarks/2ne/families/single-opposite/basis.md](sources/benchmarks/2ne/families/single-opposite/basis.md).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [benchmarks/self-energy/pi1-chain/prepared/normalization.json](sources/benchmarks/self-energy/pi1-chain/prepared/normalization.json).
Source dependency: [benchmarks/self-energy/pi1-chain/workspaces/2e-Pi1+Pi1-01/derivation.md](sources/benchmarks/self-energy/pi1-chain/workspaces/2e-Pi1+Pi1-01/derivation.md).

## 2e-Pi1+Pi1-02

Source: [benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-02/r002.json](sources/benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-02/r002.json) (revision r002); SHA-256 `60af70d4a696f29cde70d92c03775f75bd7b83ce36fd57d6583f2baee5100d60`.

Native unweighted field: `/laurent/coefficients`. unweighted scalar I; chain color multipliers are already included in the consumed channel fields.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
-1/384
```

`eps^-2` — explicit source coefficient:

```text
rho**2/256 - 1/192 + pi**2/256
```

`eps^-1` — explicit source coefficient:

```text
-3*pi**2*H_z_m1/256 - 3*H_z_m1_0_0/128 + rho**3/512 + rho**2/64 + 3*pi**2*rho/512 - 3*zeta(3)/128 - 1/48 + 23*pi**2/1536
```

`eps^0` — explicit source coefficient:

```text
-9*pi**2*H_z_0_m1/512 - 9*H_z_0_m1_0_0/256 - 3*pi**2*H_z_m1/64 + 9*H_z_m1*zeta(3)/128 - 9*pi**2*H_z_m1_0/512 - 3*H_z_m1_0_0/32 - 9*H_z_m1_0_0_0/256 + 9*pi**2*H_z_m1_m1/256 + 9*H_z_m1_m1_0_0/128 - 3*rho**4/1024 + rho**3/128 + pi**2*rho**2/1024 + 3*rho**2/64 - 9*rho*zeta(3)/256 + 3*pi**2*rho/128 - 43*zeta(3)/384 - 1/24 + 17*pi**4/5120 + 35*pi**2/768
```

### Af contribution

Consumed field: `/color_channels/1/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
5/108
```

`eps^-2`:

```text
-5*rho**2/72 - 5*pi**2/72 + 13/162
```

`eps^-1`:

```text
5*pi**2*H_z_m1/24 + 5*H_z_m1_0_0/12 - 5*rho**3/144 - 7*rho**2/27 - 5*pi**2*rho/48 - 107*pi**2/432 + 109/324 + 5*zeta(3)/12
```

`eps^0`:

```text
5*pi**2*H_z_0_m1/16 + 5*H_z_0_m1_0_0/8 - 5*H_z_m1*zeta(3)/4 + 7*pi**2*H_z_m1/9 + 5*pi**2*H_z_m1_0/16 + 14*H_z_m1_0_0/9 + 5*H_z_m1_0_0_0/8 - 5*pi**2*H_z_m1_m1/8 - 5*H_z_m1_m1_0_0/4 + 5*rho**4/96 - 7*rho**3/54 - 161*rho**2/216 - 5*pi**2*rho**2/288 - 7*pi**2*rho/18 + 5*rho*zeta(3)/8 - 235*pi**2/324 - 17*pi**4/288 + 899/1458 + 203*zeta(3)/108
```

### ff contribution

Consumed field: `/color_channels/2/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
-1/54
```

`eps^-2`:

```text
rho**2/36 - 2/81 + pi**2/36
```

`eps^-1`:

```text
-pi**2*H_z_m1/12 - H_z_m1_0_0/6 + rho**3/72 + 5*rho**2/54 + pi**2*rho/24 - zeta(3)/6 - 19/162 + 19*pi**2/216
```

`eps^0`:

```text
-pi**2*H_z_0_m1/8 - H_z_0_m1_0_0/4 - 5*pi**2*H_z_m1/18 + H_z_m1*zeta(3)/2 - pi**2*H_z_m1_0/8 - 5*H_z_m1_0_0/9 - H_z_m1_0_0_0/4 + pi**2*H_z_m1_m1/4 + H_z_m1_m1_0_0/2 - rho**4/48 + 5*rho**3/108 + pi**2*rho**2/144 + rho**2/4 - rho*zeta(3)/4 + 5*pi**2*rho/36 - 37*zeta(3)/54 - 133/729 + 17*pi**4/720 + 79*pi**2/324
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "Mkn": "existing evaluated M_{k,n}; bound Hlog files in dependencies",
    "H0": "z^(-1/2*delta)*(-z - 1)*MeijerG[{{0, 0, delta}, {}}, {{0, delta, delta}, {}}, z]/(Gamma[1 - delta]*Gamma[delta + 1]^2)"
  }
}
```

Source dependency: [benchmarks/2ne/families/single-opposite/basis-functions.json](sources/benchmarks/2ne/families/single-opposite/basis-functions.json).
Source dependency: [benchmarks/2ne/families/single-opposite/basis.md](sources/benchmarks/2ne/families/single-opposite/basis.md).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [benchmarks/self-energy/pi1-chain/prepared/normalization.json](sources/benchmarks/self-energy/pi1-chain/prepared/normalization.json).
Source dependency: [benchmarks/self-energy/pi1-chain/workspaces/2e-Pi1+Pi1-02/derivation.md](sources/benchmarks/self-energy/pi1-chain/workspaces/2e-Pi1+Pi1-02/derivation.md).

## 2e-Pi2-01

Source: [benchmarks/self-energy/pi2/results/2e-Pi2-01/r002.json](sources/benchmarks/self-energy/pi2/results/2e-Pi2-01/r002.json) (revision r002); SHA-256 `27d71a2e65cf149854fd9daa1ccc39aff7b8b06401c053d18885cf308b29b42a`.

Native unweighted field: `/I_geom_laurent/coefficients`. unweighted geometry I_geom; channel operators are already included in the consumed channel fields.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
1/576
```

`eps^-3` — explicit source coefficient:

```text
5/1152
```

`eps^-2` — explicit source coefficient:

```text
pi**2/2304 + 1/64
```

`eps^-1` — explicit source coefficient:

```text
5*pi**2/4608 + 7*zeta(3)/576 + 5/144
```

`eps^0` — explicit source coefficient:

```text
67*pi**4/276480 + 35*zeta(3)/1152 + pi**2/256 + 7/72
```

### Af contribution

Consumed field: `/color_channels/1/laurent/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
5/432
```

`eps^-3`:

```text
7/324
```

`eps^-2`:

```text
5*pi**2/1728 + zeta(3)/18 + 199/1944
```

`eps^-1`:

```text
7*pi**2/1296 + 7*zeta(3)/144 + pi**4/1080 + 4751/23328
```

`eps^0`:

```text
223*pi**4/207360 + pi**2*zeta(3)/72 + 199*pi**2/7776 + 7*zeta(5)/18 + 73*zeta(3)/216 + 97655/139968
```

### Ff contribution

Consumed field: `/color_channels/2/laurent/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
1/72
```

`eps^-2`:

```text
23/216 - zeta(3)/9
```

`eps^-1`:

```text
-5*zeta(3)/27 - pi**4/540 + pi**2/288 + 331/648
```

`eps^0`:

```text
-557*zeta(3)/648 - 7*zeta(5)/9 - pi**2*zeta(3)/36 - pi**4/324 + 23*pi**2/864 + 7993/3888
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "Mkn": "existing evaluated M_{k,n}; bound Hlog files in dependencies",
    "H0": "z^(-1/2*delta)*(-z - 1)*MeijerG[{{0, 0, delta}, {}}, {{0, delta, delta}, {}}, z]/(Gamma[1 - delta]*Gamma[delta + 1]^2)"
  }
}
```

Source dependency: [benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-01/r001.json](sources/benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-01/r001.json).
Source dependency: [benchmarks/self-energy/pi1-chain/prepared/2e-Pi1+Pi1.json](sources/benchmarks/self-energy/pi1-chain/prepared/2e-Pi1+Pi1.json).
Source dependency: [workflow/checks/self-energy-pi1-chain-solve/compute.py](sources/workflow/checks/self-energy-pi1-chain-solve/compute.py).
Source dependency: [workflow/checks/self-energy-pi1-chain-solve/exact-checks.json](sources/workflow/checks/self-energy-pi1-chain-solve/exact-checks.json).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [benchmarks/self-energy/pi2/prepared/normalization.json](sources/benchmarks/self-energy/pi2/prepared/normalization.json).
Source dependency: [benchmarks/self-energy/pi2/theory.md](sources/benchmarks/self-energy/pi2/theory.md).
Source dependency: [benchmarks/self-energy/pi2/operator-expansion-r001.json](sources/benchmarks/self-energy/pi2/operator-expansion-r001.json).
Source dependency: [benchmarks/self-energy/pi2/workspaces/2e-Pi2-01/derivation.md](sources/benchmarks/self-energy/pi2/workspaces/2e-Pi2-01/derivation.md).

A historical dependency digest differs from the current documentation file. Both digests are recorded in the JSON; no identity is inferred between the two versions. All selected coefficient-file digests match the fresh replay.


## 2e-Pi2-02

Source: [benchmarks/self-energy/pi2/results/2e-Pi2-02/r002.json](sources/benchmarks/self-energy/pi2/results/2e-Pi2-02/r002.json) (revision r002); SHA-256 `384da32d14a664a819a959f7e4abf2b0668c0e31ccfdc2963f4eaedbc9af494d`.

Native unweighted field: `/I_geom_laurent/coefficients`. unweighted geometry I_geom; channel operators are already included in the consumed channel fields.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
1/384
```

`eps^-2` — explicit source coefficient:

```text
-rho**2/256 - pi**2/256 + 1/192
```

`eps^-1` — explicit source coefficient:

```text
3*pi**2*H_z_m1/256 + 3*H_z_m1_0_0/128 - rho**3/512 - rho**2/64 - 3*pi**2*rho/512 - 23*pi**2/1536 + 1/48 + 3*zeta(3)/128
```

`eps^0` — explicit source coefficient:

```text
9*pi**2*H_z_0_m1/512 + 9*H_z_0_m1_0_0/256 - 9*H_z_m1*zeta(3)/128 + 3*pi**2*H_z_m1/64 + 9*pi**2*H_z_m1_0/512 + 3*H_z_m1_0_0/32 + 9*H_z_m1_0_0_0/256 - 9*pi**2*H_z_m1_m1/256 - 9*H_z_m1_m1_0_0/128 + 3*rho**4/1024 - rho**3/128 - 3*rho**2/64 - pi**2*rho**2/1024 - 3*pi**2*rho/128 + 9*rho*zeta(3)/256 - 35*pi**2/768 - 17*pi**4/5120 + 1/24 + 43*zeta(3)/384
```

### Af contribution

Consumed field: `/color_channels/1/laurent/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
5/288
```

`eps^-2`:

```text
-5*rho**2/192 - 5*pi**2/192 + 41/1728
```

`eps^-1`:

```text
5*pi**2*H_z_m1/64 + 5*H_z_m1_0_0/32 - 5*rho**3/384 - 101*rho**2/1152 - 5*pi**2*rho/128 - pi**2/12 + 1469/10368 + 23*zeta(3)/96
```

`eps^0`:

```text
15*pi**2*H_z_0_m1/128 + 15*H_z_0_m1_0_0/64 - 15*H_z_m1*zeta(3)/32 + 101*pi**2*H_z_m1/384 + 15*pi**2*H_z_m1_0/128 + 101*H_z_m1_0_0/192 + 15*H_z_m1_0_0_0/64 - 15*pi**2*H_z_m1_m1/64 - 15*H_z_m1_m1_0_0/32 + 5*rho**4/256 - 101*rho**3/2304 - 1961*rho**2/6912 - rho**2*zeta(3)/8 - 5*pi**2*rho**2/768 - 101*pi**2*rho/768 + 15*rho*zeta(3)/64 - 5*pi**2/18 - 239*pi**4/11520 - pi**2*zeta(3)/8 + 14597/62208 + 107*zeta(3)/192
```

### Ff contribution

Consumed field: `/color_channels/2/laurent/coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

The consumed field already includes the channel color operator; do not multiply it again.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
0
```

`eps^-2`:

```text
1/48
```

`eps^-1`:

```text
-rho**2/32 - pi**2/32 - zeta(3)/6 + 43/288
```

`eps^0`:

```text
3*pi**2*H_z_m1/32 + 3*H_z_m1_0_0/16 - rho**3/64 - 55*rho**2/192 + rho**2*zeta(3)/4 - 3*pi**2*rho/64 - 9*pi**2/32 - pi**4/360 - zeta(3)/144 + 1195/1728 + pi**2*zeta(3)/4
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "functions": {
    "L": "log(z)",
    "Mkn": "existing evaluated M_{k,n}; bound Hlog files in dependencies",
    "H0": "z^(-1/2*delta)*(-z - 1)*MeijerG[{{0, 0, delta}, {}}, {{0, delta, delta}, {}}, z]/(Gamma[1 - delta]*Gamma[delta + 1]^2)"
  }
}
```

Source dependency: [benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-02/r001.json](sources/benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-02/r001.json).
Source dependency: [benchmarks/self-energy/pi1-chain/prepared/2e-Pi1+Pi1.json](sources/benchmarks/self-energy/pi1-chain/prepared/2e-Pi1+Pi1.json).
Source dependency: [workflow/checks/self-energy-pi1-chain-solve/compute.py](sources/workflow/checks/self-energy-pi1-chain-solve/compute.py).
Source dependency: [workflow/checks/self-energy-pi1-chain-solve/exact-checks.json](sources/workflow/checks/self-energy-pi1-chain-solve/exact-checks.json).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt).
Source dependency: [workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt](sources/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt).
Source dependency: [benchmarks/self-energy/pi2/prepared/normalization.json](sources/benchmarks/self-energy/pi2/prepared/normalization.json).
Source dependency: [benchmarks/self-energy/pi2/theory.md](sources/benchmarks/self-energy/pi2/theory.md).
Source dependency: [benchmarks/self-energy/pi2/operator-expansion-r001.json](sources/benchmarks/self-energy/pi2/operator-expansion-r001.json).
Source dependency: [benchmarks/self-energy/pi2/workspaces/2e-Pi2-02/derivation.md](sources/benchmarks/self-energy/pi2/workspaces/2e-Pi2-02/derivation.md).

A historical dependency digest differs from the current documentation file. Both digests are recorded in the JSON; no identity is inferred between the two versions. All selected coefficient-file digests match the fresh replay.


## 3e-3g+Pi1-01

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-01/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-01/r003.json) (revision r003); SHA-256 `b0b007dd4d387ffd744376bb54d62dadbaaf9948421d2f3819af6d09c5200e0f`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/1152
```

`eps^-4` — explicit source coefficient:

```text
-1/1536
```

`eps^-3` — explicit source coefficient:

```text
-73/9216 + 5*pi**2/4608
```

`eps^-2` — explicit source coefficient:

```text
-527/18432 + 13*zeta(3)/1152 + 41*pi**2/18432
```

`eps^-1` — explicit source coefficient:

```text
-1447/12288 + 7*pi**2/36864 + 43*zeta(3)/1536 + 29*pi**4/36864
```

`eps^0` — explicit source coefficient:

```text
-26863/73728 + 167*zeta(3)/9216 + 91*pi**2/24576 + 47*zeta(5)/640 + 41*pi**2*zeta(3)/4608 + 4519*pi**4/2211840
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
-1/432
```

`eps^-4`:

```text
13/5184
```

`eps^-3`:

```text
655/31104 - 5*pi**2/1728
```

`eps^-2`:

```text
-103*pi**2/20736 - 13*zeta(3)/432 + 12907/186624
```

`eps^-1`:

```text
-29*pi**4/13824 - 335*zeta(3)/5184 + 263*pi**2/124416 + 317875/1119744
```

`eps^0`:

```text
-1313*pi**4/276480 - 41*pi**2*zeta(3)/1728 - 47*zeta(5)/240 - 5941*pi**2/746496 - 521*zeta(3)/31104 + 5689483/6718464
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-01/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-01/derivation.md).

## 3e-3g+Pi1-02

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-02/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-02/r003.json) (revision r003); SHA-256 `809e43e98fb22c31c63ddb4985e5c1b08d4f060a72008e28f83918daccc97389`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/2304
```

`eps^-4` — explicit source coefficient:

```text
-1/4608
```

`eps^-3` — explicit source coefficient:

```text
-65/9216 + pi**2/1024
```

`eps^-2` — explicit source coefficient:

```text
-55/2048 + 13*zeta(3)/2304 + 31*pi**2/18432
```

`eps^-1` — explicit source coefficient:

```text
-3925/36864 + 47*zeta(3)/4608 + 7*pi**2/4096 + 199*pi**4/221184
```

`eps^0` — explicit source coefficient:

```text
-25199/73728 + 17*pi**2/73728 + 127*zeta(3)/9216 + 47*zeta(5)/1280 + 23*pi**2*zeta(3)/3072 + 259*pi**4/147456
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
-1/864
```

`eps^-4`:

```text
5/5184
```

`eps^-3`:

```text
587/31104 - pi**2/384
```

`eps^-2`:

```text
-25*pi**2/6912 - 13*zeta(3)/864 + 12203/186624
```

`eps^-1`:

```text
-199*pi**4/82944 - 115*zeta(3)/5184 - 103*pi**2/41472 + 286547/1119744
```

`eps^0`:

```text
-1933*pi**4/497664 - 23*pi**2*zeta(3)/1152 - 47*zeta(5)/480 - 757*zeta(3)/31104 + 569*pi**2/248832 + 5361995/6718464
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-02/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-02/derivation.md).

## 3e-3g+Pi1-03

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-03/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-03/r003.json) (revision r003); SHA-256 `ed477d05e9da0bead80c34bb4babb08c164ca57206cae260a901def6a5e177e1`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/1152
```

`eps^-4` — explicit source coefficient:

```text
-7/2304
```

`eps^-3` — explicit source coefficient:

```text
-1/192 + 7*pi**2/4608
```

`eps^-2` — explicit source coefficient:

```text
-17/576 - pi**2/9216 + 7*zeta(3)/1152
```

`eps^-1` — explicit source coefficient:

```text
-1/18 - 49*zeta(3)/2304 + pi**2/256 + 547*pi**4/552960
```

`eps^0` — explicit source coefficient:

```text
-3/16 - 11*pi**2/2304 - 7*zeta(3)/192 + 11*pi**4/1105920 + 71*zeta(5)/1920 + 49*pi**2*zeta(3)/4608
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
-1/432
```

`eps^-4`:

```text
23/2592
```

`eps^-3`:

```text
91/7776 - 7*pi**2/1728
```

`eps^-2`:

```text
-7*zeta(3)/432 + 17*pi**2/10368 + 847/11664
```

`eps^-1`:

```text
-547*pi**4/207360 - 299*pi**2/31104 + 161*zeta(3)/2592 + 1031/8748
```

`eps^0`:

```text
-49*pi**2*zeta(3)/1728 - 71*zeta(5)/720 + 1061*pi**4/1244160 + 637*zeta(3)/7776 + 781*pi**2/46656 + 2824/6561
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-03/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-03/derivation.md).

## 3e-3g+Pi1-04

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-04/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-04/r003.json) (revision r003); SHA-256 `0f3e9e0b82b51f3821a4f612a0cefeff1eb6bff2513ad1e885fabed2a8a630ec`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/1152
```

`eps^-4` — explicit source coefficient:

```text
-7/2304
```

`eps^-3` — explicit source coefficient:

```text
-1/192 + 7*pi**2/4608
```

`eps^-2` — explicit source coefficient:

```text
-17/576 - pi**2/9216 + 7*zeta(3)/1152
```

`eps^-1` — explicit source coefficient:

```text
-1/18 - 49*zeta(3)/2304 + pi**2/256 + 547*pi**4/552960
```

`eps^0` — explicit source coefficient:

```text
-3/16 - 11*pi**2/2304 - 7*zeta(3)/192 + 11*pi**4/1105920 + 71*zeta(5)/1920 + 49*pi**2*zeta(3)/4608
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
-1/432
```

`eps^-4`:

```text
23/2592
```

`eps^-3`:

```text
91/7776 - 7*pi**2/1728
```

`eps^-2`:

```text
-7*zeta(3)/432 + 17*pi**2/10368 + 847/11664
```

`eps^-1`:

```text
-547*pi**4/207360 - 299*pi**2/31104 + 161*zeta(3)/2592 + 1031/8748
```

`eps^0`:

```text
-49*pi**2*zeta(3)/1728 - 71*zeta(5)/720 + 1061*pi**4/1244160 + 637*zeta(3)/7776 + 781*pi**2/46656 + 2824/6561
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-04/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-04/derivation.md).

## 3e-3g+Pi1-05

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-05/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-05/r003.json) (revision r003); SHA-256 `32b15a8aae68502fa406e1b659c06cfbf07cb23edbcbb0d151da9276f09c757d`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/2304
```

`eps^-4` — explicit source coefficient:

```text
-1/4608
```

`eps^-3` — explicit source coefficient:

```text
-65/9216 + pi**2/1024
```

`eps^-2` — explicit source coefficient:

```text
-55/2048 + 13*zeta(3)/2304 + 31*pi**2/18432
```

`eps^-1` — explicit source coefficient:

```text
-3925/36864 + 47*zeta(3)/4608 + 7*pi**2/4096 + 199*pi**4/221184
```

`eps^0` — explicit source coefficient:

```text
-25199/73728 + 17*pi**2/73728 + 127*zeta(3)/9216 + 47*zeta(5)/1280 + 23*pi**2*zeta(3)/3072 + 259*pi**4/147456
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
-1/864
```

`eps^-4`:

```text
5/5184
```

`eps^-3`:

```text
587/31104 - pi**2/384
```

`eps^-2`:

```text
-25*pi**2/6912 - 13*zeta(3)/864 + 12203/186624
```

`eps^-1`:

```text
-199*pi**4/82944 - 115*zeta(3)/5184 - 103*pi**2/41472 + 286547/1119744
```

`eps^0`:

```text
-1933*pi**4/497664 - 23*pi**2*zeta(3)/1152 - 47*zeta(5)/480 - 757*zeta(3)/31104 + 569*pi**2/248832 + 5361995/6718464
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-05/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-05/derivation.md).

## 3e-3g+Pi1-06

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-06/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-06/r003.json) (revision r003); SHA-256 `8a334297e7e8afe8b95f0bed3c24fc826ec67b13792df8f29105b27cc0dbfcd5`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/1152
```

`eps^-4` — explicit source coefficient:

```text
-1/1536
```

`eps^-3` — explicit source coefficient:

```text
-73/9216 + 5*pi**2/4608
```

`eps^-2` — explicit source coefficient:

```text
-527/18432 + 13*zeta(3)/1152 + 41*pi**2/18432
```

`eps^-1` — explicit source coefficient:

```text
-1447/12288 + 7*pi**2/36864 + 43*zeta(3)/1536 + 29*pi**4/36864
```

`eps^0` — explicit source coefficient:

```text
-26863/73728 + 167*zeta(3)/9216 + 91*pi**2/24576 + 47*zeta(5)/640 + 41*pi**2*zeta(3)/4608 + 4519*pi**4/2211840
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
-1/432
```

`eps^-4`:

```text
13/5184
```

`eps^-3`:

```text
655/31104 - 5*pi**2/1728
```

`eps^-2`:

```text
-103*pi**2/20736 - 13*zeta(3)/432 + 12907/186624
```

`eps^-1`:

```text
-29*pi**4/13824 - 335*zeta(3)/5184 + 263*pi**2/124416 + 317875/1119744
```

`eps^0`:

```text
-1313*pi**4/276480 - 41*pi**2*zeta(3)/1728 - 47*zeta(5)/240 - 5941*pi**2/746496 - 521*zeta(3)/31104 + 5689483/6718464
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-06/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-06/derivation.md).

## 3e-3g+Pi1-07

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-07/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-07/r003.json) (revision r003); SHA-256 `fac4e545879c642b4c0eb3214760e9ba467ecc14ac522f1ee842bc3b7389dcb0`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — zero established by source pole metadata:

```text
0
```

`eps^-4` — zero established by source pole metadata:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
1/768
```

`eps^-2` — explicit source coefficient:

```text
pi**2*H_r_0/256 + H_r_0_0/64 + H_r_0_0_0/128 - H_r_0_0_1/128 - H_r_0_1/64 - H_r_0_1_0/128 + H_r_0_1_1/128 + pi**2*H_r_1/256 - H_r_1_0/64 + H_r_1_0_0/128 - H_r_1_0_1/128 + H_r_1_1/64 - H_r_1_1_0/128 + H_r_1_1_1/128 - zeta(3)/64 - 5/768 + pi**2/128
```

`eps^-1` — explicit source coefficient:

```text
-H_r_0*(pi**2*r/3072 - 65*pi**2/3072 + 3*zeta(3)/128) - H_r_0_0*(-87/1024 - pi**2/512) - H_r_0_0_0*(r/512 - 11/256) + 3*H_r_0_0_0_0/256 - 3*H_r_0_0_0_1/256 - H_r_0_0_1*(11/256 - r/512) - 3*H_r_0_0_1_0/256 + 3*H_r_0_0_1_1/256 - H_r_0_1*(87/1024 - 5*pi**2/512) - H_r_0_1_0*(11/256 - r/512) + 3*H_r_0_1_0_0/256 - 3*H_r_0_1_0_1/256 - H_r_0_1_1*(r/512 - 11/256) - 3*H_r_0_1_1_0/256 + 3*H_r_0_1_1_1/256 - H_r_1*(-pi**2*r/3072 - pi**2/48 + 3*zeta(3)/128) - H_r_1_0*(87/1024 - 5*pi**2/512) - H_r_1_0_0*(-r/512 - 21/512) + 3*H_r_1_0_0_0/256 - 3*H_r_1_0_0_1/256 - H_r_1_0_1*(r/512 + 21/512) - 3*H_r_1_0_1_0/256 + 3*H_r_1_0_1_1/256 - H_r_1_1*(-87/1024 - pi**2/512) - H_r_1_1_0*(r/512 + 21/512) + 3*H_r_1_1_0_0/256 - 3*H_r_1_1_0_1/256 - H_r_1_1_1*(-r/512 - 21/512) - 3*H_r_1_1_1_0/256 + 3*H_r_1_1_1_1/256 - 43*zeta(3)/512 - 1/192 + pi**4/10240 + 275*pi**2/6144
```

`eps^0` — explicit source coefficient:

```text
-H_r_0*(pi**2*r/3072 - 995*pi**2/12288 - 209*pi**4/61440 + 129*zeta(3)/1024) - H_r_0_0*(-667/2048 - 5*pi**2/512 + 9*zeta(3)/256) - H_r_0_0_0*(r/512 - 333/2048 + 47*pi**2/3072) - 27*H_r_0_0_0_0/256 - 71*H_r_0_0_0_0_0/512 + 71*H_r_0_0_0_0_1/512 + 27*H_r_0_0_0_1/256 + 71*H_r_0_0_0_1_0/512 - 71*H_r_0_0_0_1_1/512 - H_r_0_0_1*(-r/512 - 47*pi**2/3072 + 333/2048) + 27*H_r_0_0_1_0/256 + 71*H_r_0_0_1_0_0/512 - 71*H_r_0_0_1_0_1/512 - 27*H_r_0_0_1_1/256 - 71*H_r_0_0_1_1_0/512 + 71*H_r_0_0_1_1_1/512 - H_r_0_1*(-109*pi**2/2048 + 9*zeta(3)/256 + 667/2048) - H_r_0_1_0*(-r/512 - 47*pi**2/3072 + 333/2048) + 237*H_r_0_1_0_0/1024 + 35*H_r_0_1_0_0_0/512 - 35*H_r_0_1_0_0_1/512 - 237*H_r_0_1_0_1/1024 - 35*H_r_0_1_0_1_0/512 + 35*H_r_0_1_0_1_1/512 - H_r_0_1_1*(r/512 - 61*pi**2/3072 - 333/2048) - 237*H_r_0_1_1_0/1024 + H_r_0_1_1_0_0/512 - H_r_0_1_1_0_1/512 + 237*H_r_0_1_1_1/1024 - H_r_0_1_1_1_0/512 + H_r_0_1_1_1_1/512 - H_r_1*(-pi**2*r/3072 - 991*pi**2/12288 - 209*pi**4/61440 + 129*zeta(3)/1024) - H_r_1_0*(-109*pi**2/2048 + 9*zeta(3)/256 + 667/2048) - H_r_1_0_0*(-r/512 - 61*pi**2/3072 - 329/2048) + 237*H_r_1_0_0_0/1024 + H_r_1_0_0_0_0/512 - H_r_1_0_0_0_1/512 - 237*H_r_1_0_0_1/1024 - H_r_1_0_0_1_0/512 + H_r_1_0_0_1_1/512 - H_r_1_0_1*(r/512 - 47*pi**2/3072 + 329/2048) - 237*H_r_1_0_1_0/1024 + 35*H_r_1_0_1_0_0/512 - 35*H_r_1_0_1_0_1/512 + 237*H_r_1_0_1_1/1024 - 35*H_r_1_0_1_1_0/512 + 35*H_r_1_0_1_1_1/512 - H_r_1_1*(-667/2048 - 5*pi**2/512 + 9*zeta(3)/256) - H_r_1_1_0*(r/512 - 47*pi**2/3072 + 329/2048) - 27*H_r_1_1_0_0/256 + 71*H_r_1_1_0_0_0/512 - 71*H_r_1_1_0_0_1/512 + 27*H_r_1_1_0_1/256 - 71*H_r_1_1_0_1_0/512 + 71*H_r_1_1_0_1_1/512 - H_r_1_1_1*(-r/512 - 329/2048 + 47*pi**2/3072) + 27*H_r_1_1_1_0/256 - 71*H_r_1_1_1_0_0/512 + 71*H_r_1_1_1_0_1/512 - 27*H_r_1_1_1_1/256 + 71*H_r_1_1_1_1_0/512 - 71*H_r_1_1_1_1_1/512 - 917*zeta(3)/3072 - 31*pi**2*zeta(3)/1536 - 37*zeta(5)/256 - 35/512 + 147*pi**4/20480 + 1957*pi**2/12288
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
-1/288
```

`eps^-2`:

```text
-pi**2*H_r_0/96 - H_r_0_0/24 - H_r_0_0_0/48 + H_r_0_0_1/48 + H_r_0_1/24 + H_r_0_1_0/48 - H_r_0_1_1/48 - pi**2*H_r_1/96 + H_r_1_0/24 - H_r_1_0_0/48 + H_r_1_0_1/48 - H_r_1_1/24 + H_r_1_1_0/48 - H_r_1_1_1/48 - pi**2/48 + 1/54 + zeta(3)/24
```

`eps^-1`:

```text
pi**2*H_r_0*r/1152 - 61*pi**2*H_r_0/1152 + H_r_0*zeta(3)/16 - 245*H_r_0_0/1152 - pi**2*H_r_0_0/192 + H_r_0_0_0*r/192 - 31*H_r_0_0_0/288 - H_r_0_0_0_0/32 + H_r_0_0_0_1/32 - H_r_0_0_1*r/192 + 31*H_r_0_0_1/288 + H_r_0_0_1_0/32 - H_r_0_0_1_1/32 - 5*pi**2*H_r_0_1/192 + 245*H_r_0_1/1152 - H_r_0_1_0*r/192 + 31*H_r_0_1_0/288 - H_r_0_1_0_0/32 + H_r_0_1_0_1/32 + H_r_0_1_1*r/192 - 31*H_r_0_1_1/288 + H_r_0_1_1_0/32 - H_r_0_1_1_1/32 - pi**2*H_r_1*r/1152 - 5*pi**2*H_r_1/96 + H_r_1*zeta(3)/16 - 5*pi**2*H_r_1_0/192 + 245*H_r_1_0/1152 - H_r_1_0_0*r/192 - 59*H_r_1_0_0/576 - H_r_1_0_0_0/32 + H_r_1_0_0_1/32 + H_r_1_0_1*r/192 + 59*H_r_1_0_1/576 + H_r_1_0_1_0/32 - H_r_1_0_1_1/32 - 245*H_r_1_1/1152 - pi**2*H_r_1_1/192 + H_r_1_1_0*r/192 + 59*H_r_1_1_0/576 - H_r_1_1_0_0/32 + H_r_1_1_0_1/32 - H_r_1_1_1*r/192 - 59*H_r_1_1_1/576 + H_r_1_1_1_0/32 - H_r_1_1_1_1/32 - 259*pi**2/2304 - pi**4/3840 + 23/2592 + 121*zeta(3)/576
```

`eps^0`:

```text
pi**2*H_r_0*r/1728 - 2693*pi**2*H_r_0/13824 - 209*pi**4*H_r_0/23040 + 121*H_r_0*zeta(3)/384 - 5417*H_r_0_0/6912 - 7*pi**2*H_r_0_0/288 + 3*H_r_0_0*zeta(3)/32 + H_r_0_0_0*r/288 - 2701*H_r_0_0_0/6912 + 47*pi**2*H_r_0_0_0/1152 + 7*H_r_0_0_0_0/24 + 71*H_r_0_0_0_0_0/192 - 71*H_r_0_0_0_0_1/192 - 7*H_r_0_0_0_1/24 - 71*H_r_0_0_0_1_0/192 + 71*H_r_0_0_0_1_1/192 - H_r_0_0_1*r/288 - 47*pi**2*H_r_0_0_1/1152 + 2701*H_r_0_0_1/6912 - 7*H_r_0_0_1_0/24 - 71*H_r_0_0_1_0_0/192 + 71*H_r_0_0_1_0_1/192 + 7*H_r_0_0_1_1/24 + 71*H_r_0_0_1_1_0/192 - 71*H_r_0_0_1_1_1/192 - 307*pi**2*H_r_0_1/2304 + 3*H_r_0_1*zeta(3)/32 + 5417*H_r_0_1/6912 - H_r_0_1_0*r/288 - 47*pi**2*H_r_0_1_0/1152 + 2701*H_r_0_1_0/6912 - 233*H_r_0_1_0_0/384 - 35*H_r_0_1_0_0_0/192 + 35*H_r_0_1_0_0_1/192 + 233*H_r_0_1_0_1/384 + 35*H_r_0_1_0_1_0/192 - 35*H_r_0_1_0_1_1/192 + H_r_0_1_1*r/288 - 61*pi**2*H_r_0_1_1/1152 - 2701*H_r_0_1_1/6912 + 233*H_r_0_1_1_0/384 - H_r_0_1_1_0_0/192 + H_r_0_1_1_0_1/192 - 233*H_r_0_1_1_1/384 + H_r_0_1_1_1_0/192 - H_r_0_1_1_1_1/192 - pi**2*H_r_1*r/1728 - 895*pi**2*H_r_1/4608 - 209*pi**4*H_r_1/23040 + 121*H_r_1*zeta(3)/384 - 307*pi**2*H_r_1_0/2304 + 3*H_r_1_0*zeta(3)/32 + 5417*H_r_1_0/6912 - H_r_1_0_0*r/288 - 61*pi**2*H_r_1_0_0/1152 - 2677*H_r_1_0_0/6912 - 233*H_r_1_0_0_0/384 - H_r_1_0_0_0_0/192 + H_r_1_0_0_0_1/192 + 233*H_r_1_0_0_1/384 + H_r_1_0_0_1_0/192 - H_r_1_0_0_1_1/192 + H_r_1_0_1*r/288 - 47*pi**2*H_r_1_0_1/1152 + 2677*H_r_1_0_1/6912 + 233*H_r_1_0_1_0/384 - 35*H_r_1_0_1_0_0/192 + 35*H_r_1_0_1_0_1/192 - 233*H_r_1_0_1_1/384 + 35*H_r_1_0_1_1_0/192 - 35*H_r_1_0_1_1_1/192 - 5417*H_r_1_1/6912 - 7*pi**2*H_r_1_1/288 + 3*H_r_1_1*zeta(3)/32 + H_r_1_1_0*r/288 - 47*pi**2*H_r_1_1_0/1152 + 2677*H_r_1_1_0/6912 + 7*H_r_1_1_0_0/24 - 71*H_r_1_1_0_0_0/192 + 71*H_r_1_1_0_0_1/192 - 7*H_r_1_1_0_1/24 + 71*H_r_1_1_0_1_0/192 - 71*H_r_1_1_0_1_1/192 - H_r_1_1_1*r/288 - 2677*H_r_1_1_1/6912 + 47*pi**2*H_r_1_1_1/1152 - 7*H_r_1_1_1_0/24 + 71*H_r_1_1_1_0_0/192 - 71*H_r_1_1_1_0_1/192 + 7*H_r_1_1_1_1/24 - 71*H_r_1_1_1_1_0/192 + 71*H_r_1_1_1_1_1/192 - 5257*pi**2/13824 - 439*pi**4/23040 + 2711/15552 + 37*zeta(5)/96 + 31*pi**2*zeta(3)/576 + 2461*zeta(3)/3456
```

Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "expression_language": "Maple/HyperInt; zeta[n]=Zeta(n); Hlog uses dt/(t-letter)",
  "variables": {
    "r": "1/(1+z)",
    "z": "positive real kinematic ratio; z1 in the raw source"
  },
  "Hlog_convention": {
    "empty_word": "1",
    "recursion": "Hlog(x,[a,w])=integral_0^x dt/(t-a) Hlog(t,w)",
    "all_zero_word": "log(x)^n/n!",
    "Hlog(x,[1])": "log(1-x)",
    "branch": "real tangential-basepoint GPL on 0<x<1"
  },
  "M_definition": "M(k,j,z)=integral_0^1 [log(beta*(1-beta))]^k [log(z)+log((1-beta)/beta)]^j/(beta-z/(1+z)) d beta, j>=1; removable internal point"
}
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-07/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-07/derivation.md).

## 3e-3g+Pi1-08

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-08/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-08/r003.json) (revision r003); SHA-256 `cb3f9d781bea07cb85fb9547a559dd6c59f8f044c9fe2ceb7b5434b9ab75d801`.

Native unweighted field: `/laurent_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — zero established by source pole metadata:

```text
0
```

`eps^-4` — zero established by source pole metadata:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
1/768
```

`eps^-2` — explicit source coefficient:

```text
pi**2*H_r_0/256 + H_r_0_0/64 + H_r_0_0_0/128 - H_r_0_0_1/128 - H_r_0_1/64 - H_r_0_1_0/128 + H_r_0_1_1/128 + pi**2*H_r_1/256 - H_r_1_0/64 + H_r_1_0_0/128 - H_r_1_0_1/128 + H_r_1_1/64 - H_r_1_1_0/128 + H_r_1_1_1/128 - zeta(3)/64 - 5/768 + pi**2/128
```

`eps^-1` — explicit source coefficient:

```text
-H_r_0*(pi**2*r/3072 - 65*pi**2/3072 + 3*zeta(3)/128) - H_r_0_0*(-87/1024 - pi**2/512) - H_r_0_0_0*(r/512 - 11/256) + 3*H_r_0_0_0_0/256 - 3*H_r_0_0_0_1/256 - H_r_0_0_1*(11/256 - r/512) - 3*H_r_0_0_1_0/256 + 3*H_r_0_0_1_1/256 - H_r_0_1*(87/1024 - 5*pi**2/512) - H_r_0_1_0*(11/256 - r/512) + 3*H_r_0_1_0_0/256 - 3*H_r_0_1_0_1/256 - H_r_0_1_1*(r/512 - 11/256) - 3*H_r_0_1_1_0/256 + 3*H_r_0_1_1_1/256 - H_r_1*(-pi**2*r/3072 - pi**2/48 + 3*zeta(3)/128) - H_r_1_0*(87/1024 - 5*pi**2/512) - H_r_1_0_0*(-r/512 - 21/512) + 3*H_r_1_0_0_0/256 - 3*H_r_1_0_0_1/256 - H_r_1_0_1*(r/512 + 21/512) - 3*H_r_1_0_1_0/256 + 3*H_r_1_0_1_1/256 - H_r_1_1*(-87/1024 - pi**2/512) - H_r_1_1_0*(r/512 + 21/512) + 3*H_r_1_1_0_0/256 - 3*H_r_1_1_0_1/256 - H_r_1_1_1*(-r/512 - 21/512) - 3*H_r_1_1_1_0/256 + 3*H_r_1_1_1_1/256 - 43*zeta(3)/512 - 1/192 + pi**4/10240 + 275*pi**2/6144
```

`eps^0` — explicit source coefficient:

```text
-H_r_0*(pi**2*r/3072 - 995*pi**2/12288 - 209*pi**4/61440 + 129*zeta(3)/1024) - H_r_0_0*(-667/2048 - 5*pi**2/512 + 9*zeta(3)/256) - H_r_0_0_0*(r/512 - 333/2048 + 47*pi**2/3072) - 27*H_r_0_0_0_0/256 - 71*H_r_0_0_0_0_0/512 + 71*H_r_0_0_0_0_1/512 + 27*H_r_0_0_0_1/256 + 71*H_r_0_0_0_1_0/512 - 71*H_r_0_0_0_1_1/512 - H_r_0_0_1*(-r/512 - 47*pi**2/3072 + 333/2048) + 27*H_r_0_0_1_0/256 + 71*H_r_0_0_1_0_0/512 - 71*H_r_0_0_1_0_1/512 - 27*H_r_0_0_1_1/256 - 71*H_r_0_0_1_1_0/512 + 71*H_r_0_0_1_1_1/512 - H_r_0_1*(-109*pi**2/2048 + 9*zeta(3)/256 + 667/2048) - H_r_0_1_0*(-r/512 - 47*pi**2/3072 + 333/2048) + 237*H_r_0_1_0_0/1024 + 35*H_r_0_1_0_0_0/512 - 35*H_r_0_1_0_0_1/512 - 237*H_r_0_1_0_1/1024 - 35*H_r_0_1_0_1_0/512 + 35*H_r_0_1_0_1_1/512 - H_r_0_1_1*(r/512 - 61*pi**2/3072 - 333/2048) - 237*H_r_0_1_1_0/1024 + H_r_0_1_1_0_0/512 - H_r_0_1_1_0_1/512 + 237*H_r_0_1_1_1/1024 - H_r_0_1_1_1_0/512 + H_r_0_1_1_1_1/512 - H_r_1*(-pi**2*r/3072 - 991*pi**2/12288 - 209*pi**4/61440 + 129*zeta(3)/1024) - H_r_1_0*(-109*pi**2/2048 + 9*zeta(3)/256 + 667/2048) - H_r_1_0_0*(-r/512 - 61*pi**2/3072 - 329/2048) + 237*H_r_1_0_0_0/1024 + H_r_1_0_0_0_0/512 - H_r_1_0_0_0_1/512 - 237*H_r_1_0_0_1/1024 - H_r_1_0_0_1_0/512 + H_r_1_0_0_1_1/512 - H_r_1_0_1*(r/512 - 47*pi**2/3072 + 329/2048) - 237*H_r_1_0_1_0/1024 + 35*H_r_1_0_1_0_0/512 - 35*H_r_1_0_1_0_1/512 + 237*H_r_1_0_1_1/1024 - 35*H_r_1_0_1_1_0/512 + 35*H_r_1_0_1_1_1/512 - H_r_1_1*(-667/2048 - 5*pi**2/512 + 9*zeta(3)/256) - H_r_1_1_0*(r/512 - 47*pi**2/3072 + 329/2048) - 27*H_r_1_1_0_0/256 + 71*H_r_1_1_0_0_0/512 - 71*H_r_1_1_0_0_1/512 + 27*H_r_1_1_0_1/256 - 71*H_r_1_1_0_1_0/512 + 71*H_r_1_1_0_1_1/512 - H_r_1_1_1*(-r/512 - 329/2048 + 47*pi**2/3072) + 27*H_r_1_1_1_0/256 - 71*H_r_1_1_1_0_0/512 + 71*H_r_1_1_1_0_1/512 - 27*H_r_1_1_1_1/256 + 71*H_r_1_1_1_1_0/512 - 71*H_r_1_1_1_1_1/512 - 917*zeta(3)/3072 - 31*pi**2*zeta(3)/1536 - 37*zeta(5)/256 - 35/512 + 147*pi**4/20480 + 1957*pi**2/12288
```

### Af contribution

Consumed field: `/laurent_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
-1/288
```

`eps^-2`:

```text
-pi**2*H_r_0/96 - H_r_0_0/24 - H_r_0_0_0/48 + H_r_0_0_1/48 + H_r_0_1/24 + H_r_0_1_0/48 - H_r_0_1_1/48 - pi**2*H_r_1/96 + H_r_1_0/24 - H_r_1_0_0/48 + H_r_1_0_1/48 - H_r_1_1/24 + H_r_1_1_0/48 - H_r_1_1_1/48 - pi**2/48 + 1/54 + zeta(3)/24
```

`eps^-1`:

```text
pi**2*H_r_0*r/1152 - 61*pi**2*H_r_0/1152 + H_r_0*zeta(3)/16 - 245*H_r_0_0/1152 - pi**2*H_r_0_0/192 + H_r_0_0_0*r/192 - 31*H_r_0_0_0/288 - H_r_0_0_0_0/32 + H_r_0_0_0_1/32 - H_r_0_0_1*r/192 + 31*H_r_0_0_1/288 + H_r_0_0_1_0/32 - H_r_0_0_1_1/32 - 5*pi**2*H_r_0_1/192 + 245*H_r_0_1/1152 - H_r_0_1_0*r/192 + 31*H_r_0_1_0/288 - H_r_0_1_0_0/32 + H_r_0_1_0_1/32 + H_r_0_1_1*r/192 - 31*H_r_0_1_1/288 + H_r_0_1_1_0/32 - H_r_0_1_1_1/32 - pi**2*H_r_1*r/1152 - 5*pi**2*H_r_1/96 + H_r_1*zeta(3)/16 - 5*pi**2*H_r_1_0/192 + 245*H_r_1_0/1152 - H_r_1_0_0*r/192 - 59*H_r_1_0_0/576 - H_r_1_0_0_0/32 + H_r_1_0_0_1/32 + H_r_1_0_1*r/192 + 59*H_r_1_0_1/576 + H_r_1_0_1_0/32 - H_r_1_0_1_1/32 - 245*H_r_1_1/1152 - pi**2*H_r_1_1/192 + H_r_1_1_0*r/192 + 59*H_r_1_1_0/576 - H_r_1_1_0_0/32 + H_r_1_1_0_1/32 - H_r_1_1_1*r/192 - 59*H_r_1_1_1/576 + H_r_1_1_1_0/32 - H_r_1_1_1_1/32 - 259*pi**2/2304 - pi**4/3840 + 23/2592 + 121*zeta(3)/576
```

`eps^0`:

```text
pi**2*H_r_0*r/1728 - 2693*pi**2*H_r_0/13824 - 209*pi**4*H_r_0/23040 + 121*H_r_0*zeta(3)/384 - 5417*H_r_0_0/6912 - 7*pi**2*H_r_0_0/288 + 3*H_r_0_0*zeta(3)/32 + H_r_0_0_0*r/288 - 2701*H_r_0_0_0/6912 + 47*pi**2*H_r_0_0_0/1152 + 7*H_r_0_0_0_0/24 + 71*H_r_0_0_0_0_0/192 - 71*H_r_0_0_0_0_1/192 - 7*H_r_0_0_0_1/24 - 71*H_r_0_0_0_1_0/192 + 71*H_r_0_0_0_1_1/192 - H_r_0_0_1*r/288 - 47*pi**2*H_r_0_0_1/1152 + 2701*H_r_0_0_1/6912 - 7*H_r_0_0_1_0/24 - 71*H_r_0_0_1_0_0/192 + 71*H_r_0_0_1_0_1/192 + 7*H_r_0_0_1_1/24 + 71*H_r_0_0_1_1_0/192 - 71*H_r_0_0_1_1_1/192 - 307*pi**2*H_r_0_1/2304 + 3*H_r_0_1*zeta(3)/32 + 5417*H_r_0_1/6912 - H_r_0_1_0*r/288 - 47*pi**2*H_r_0_1_0/1152 + 2701*H_r_0_1_0/6912 - 233*H_r_0_1_0_0/384 - 35*H_r_0_1_0_0_0/192 + 35*H_r_0_1_0_0_1/192 + 233*H_r_0_1_0_1/384 + 35*H_r_0_1_0_1_0/192 - 35*H_r_0_1_0_1_1/192 + H_r_0_1_1*r/288 - 61*pi**2*H_r_0_1_1/1152 - 2701*H_r_0_1_1/6912 + 233*H_r_0_1_1_0/384 - H_r_0_1_1_0_0/192 + H_r_0_1_1_0_1/192 - 233*H_r_0_1_1_1/384 + H_r_0_1_1_1_0/192 - H_r_0_1_1_1_1/192 - pi**2*H_r_1*r/1728 - 895*pi**2*H_r_1/4608 - 209*pi**4*H_r_1/23040 + 121*H_r_1*zeta(3)/384 - 307*pi**2*H_r_1_0/2304 + 3*H_r_1_0*zeta(3)/32 + 5417*H_r_1_0/6912 - H_r_1_0_0*r/288 - 61*pi**2*H_r_1_0_0/1152 - 2677*H_r_1_0_0/6912 - 233*H_r_1_0_0_0/384 - H_r_1_0_0_0_0/192 + H_r_1_0_0_0_1/192 + 233*H_r_1_0_0_1/384 + H_r_1_0_0_1_0/192 - H_r_1_0_0_1_1/192 + H_r_1_0_1*r/288 - 47*pi**2*H_r_1_0_1/1152 + 2677*H_r_1_0_1/6912 + 233*H_r_1_0_1_0/384 - 35*H_r_1_0_1_0_0/192 + 35*H_r_1_0_1_0_1/192 - 233*H_r_1_0_1_1/384 + 35*H_r_1_0_1_1_0/192 - 35*H_r_1_0_1_1_1/192 - 5417*H_r_1_1/6912 - 7*pi**2*H_r_1_1/288 + 3*H_r_1_1*zeta(3)/32 + H_r_1_1_0*r/288 - 47*pi**2*H_r_1_1_0/1152 + 2677*H_r_1_1_0/6912 + 7*H_r_1_1_0_0/24 - 71*H_r_1_1_0_0_0/192 + 71*H_r_1_1_0_0_1/192 - 7*H_r_1_1_0_1/24 + 71*H_r_1_1_0_1_0/192 - 71*H_r_1_1_0_1_1/192 - H_r_1_1_1*r/288 - 2677*H_r_1_1_1/6912 + 47*pi**2*H_r_1_1_1/1152 - 7*H_r_1_1_1_0/24 + 71*H_r_1_1_1_0_0/192 - 71*H_r_1_1_1_0_1/192 + 7*H_r_1_1_1_1/24 - 71*H_r_1_1_1_1_0/192 + 71*H_r_1_1_1_1_1/192 - 5257*pi**2/13824 - 439*pi**4/23040 + 2711/15552 + 37*zeta(5)/96 + 31*pi**2*zeta(3)/576 + 2461*zeta(3)/3456
```

Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "expression_language": "Maple/HyperInt; zeta[n]=Zeta(n); Hlog uses dt/(t-letter)",
  "variables": {
    "r": "1/(1+z)",
    "z": "positive real kinematic ratio; z1 in the raw source"
  },
  "Hlog_convention": {
    "empty_word": "1",
    "recursion": "Hlog(x,[a,w])=integral_0^x dt/(t-a) Hlog(t,w)",
    "all_zero_word": "log(x)^n/n!",
    "Hlog(x,[1])": "log(1-x)",
    "branch": "real tangential-basepoint GPL on 0<x<1"
  },
  "M_definition": "M(k,j,z)=integral_0^1 [log(beta*(1-beta))]^k [log(z)+log((1-beta)/beta)]^j/(beta-z/(1+z)) d beta, j>=1; removable internal point"
}
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-08/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-08/derivation.md).

## 3e-3g+Pi1-09

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-09/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-09/r003.json) (revision r003); SHA-256 `dfc0663ca2e01e55574ddf3ae50e7149f4fd294943ffe3e842237e8ca079e326`.

Native unweighted field: `/laurent_Hlog_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — zero established by source pole metadata:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
1/384
```

`eps^-3` — explicit source coefficient:

```text
-1/96
```

`eps^-2` — explicit source coefficient:

```text
pi**2*H_r_0/128 + H_r_0_0/32 + H_r_0_0_0/64 - H_r_0_0_1/64 - H_r_0_1/32 - H_r_0_1_0/64 + H_r_0_1_1/64 + pi**2*H_r_1/128 - H_r_1_0/32 + H_r_1_0_0/64 - H_r_1_0_1/64 + H_r_1_1/32 - H_r_1_1_0/64 + H_r_1_1_1/64 - zeta(3)/32 - 1/96 + 31*pi**2/1536
```

`eps^-1` — explicit source coefficient:

```text
-H_r_0*(-5*pi**2/128 + zeta(3)/16) - H_r_0_0*(-1/8 - pi**2/128) + 5*H_r_0_0_0/64 + H_r_0_0_0_0/32 - H_r_0_0_0_1/32 - 5*H_r_0_0_1/64 - H_r_0_0_1_0/32 + H_r_0_0_1_1/32 - H_r_0_1*(1/8 - 3*pi**2/128) - 5*H_r_0_1_0/64 + H_r_0_1_0_0/32 - H_r_0_1_0_1/32 + 5*H_r_0_1_1/64 - H_r_0_1_1_0/32 + H_r_0_1_1_1/32 - H_r_1*(-5*pi**2/128 + zeta(3)/16) - H_r_1_0*(1/8 - 3*pi**2/128) + 5*H_r_1_0_0/64 + H_r_1_0_0_0/32 - H_r_1_0_0_1/32 - 5*H_r_1_0_1/64 - H_r_1_0_1_0/32 + H_r_1_0_1_1/32 - H_r_1_1*(-1/8 - pi**2/128) - 5*H_r_1_1_0/64 + H_r_1_1_0_0/32 - H_r_1_1_0_1/32 + 5*H_r_1_1_1/64 - H_r_1_1_1_0/32 + H_r_1_1_1_1/32 - 53*zeta(3)/384 - 1/12 + pi**4/640 + 23*pi**2/384
```

`eps^0` — explicit source coefficient:

```text
-H_r_0*(-pi**2/8 - 19*pi**4/2560 + 17*zeta(3)/64) - H_r_0_0*(-3/8 - 3*pi**2/128 + zeta(3)/32) - H_r_0_0_0*(-1/4 + 7*pi**2/256) - 7*H_r_0_0_0_0/32 - 17*H_r_0_0_0_0_0/64 + 17*H_r_0_0_0_0_1/64 + 7*H_r_0_0_0_1/32 + 17*H_r_0_0_0_1_0/64 - 17*H_r_0_0_0_1_1/64 - H_r_0_0_1*(1/4 - 11*pi**2/256) + 7*H_r_0_0_1_0/32 + 19*H_r_0_0_1_0_0/64 - 19*H_r_0_0_1_0_1/64 - 7*H_r_0_0_1_1/32 - 19*H_r_0_0_1_1_0/64 + 19*H_r_0_0_1_1_1/64 - H_r_0_1*(-7*pi**2/64 + 11*zeta(3)/64 + 3/8) - H_r_0_1_0*(1/4 - 11*pi**2/256) + 31*H_r_0_1_0_0/64 + 5*H_r_0_1_0_0_0/32 - 5*H_r_0_1_0_0_1/32 - 31*H_r_0_1_0_1/64 - 5*H_r_0_1_0_1_0/32 + 5*H_r_0_1_0_1_1/32 - H_r_0_1_1*(-11*pi**2/256 - 1/4) - 31*H_r_0_1_1_0/64 + H_r_0_1_1_0_0/64 - H_r_0_1_1_0_1/64 + 31*H_r_0_1_1_1/64 - H_r_0_1_1_1_0/64 + H_r_0_1_1_1_1/64 - H_r_1*(-pi**2/8 - 19*pi**4/2560 + 17*zeta(3)/64) - H_r_1_0*(-7*pi**2/64 + 11*zeta(3)/64 + 3/8) - H_r_1_0_0*(-11*pi**2/256 - 1/4) + 31*H_r_1_0_0_0/64 + H_r_1_0_0_0_0/64 - H_r_1_0_0_0_1/64 - 31*H_r_1_0_0_1/64 - H_r_1_0_0_1_0/64 + H_r_1_0_0_1_1/64 - H_r_1_0_1*(1/4 - 11*pi**2/256) - 31*H_r_1_0_1_0/64 + 5*H_r_1_0_1_0_0/32 - 5*H_r_1_0_1_0_1/32 + 31*H_r_1_0_1_1/64 - 5*H_r_1_0_1_1_0/32 + 5*H_r_1_0_1_1_1/32 - H_r_1_1*(-3/8 - 3*pi**2/128 + zeta(3)/32) - H_r_1_1_0*(1/4 - 11*pi**2/256) - 7*H_r_1_1_0_0/32 + 19*H_r_1_1_0_0_0/64 - 19*H_r_1_1_0_0_1/64 + 7*H_r_1_1_0_1/32 - 19*H_r_1_1_0_1_0/64 + 19*H_r_1_1_0_1_1/64 - H_r_1_1_1*(-1/4 + 7*pi**2/256) + 7*H_r_1_1_1_0/32 - 17*H_r_1_1_1_0_0/64 + 17*H_r_1_1_1_0_1/64 - 7*H_r_1_1_1_1/32 + 17*H_r_1_1_1_1_0/64 - 17*H_r_1_1_1_1_1/64 - 55*zeta(3)/96 - 3*zeta(5)/8 - pi**2*zeta(3)/64 - 1/8 + 3571*pi**4/184320 + 77*pi**2/384
```

### Af contribution

Consumed field: `/laurent_Hlog_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
-1/144
```

`eps^-3`:

```text
13/432
```

`eps^-2`:

```text
-pi**2*H_r_0/48 - H_r_0_0/12 - H_r_0_0_0/24 + H_r_0_0_1/24 + H_r_0_1/12 + H_r_0_1_0/24 - H_r_0_1_1/24 - pi**2*H_r_1/48 + H_r_1_0/12 - H_r_1_0_0/24 + H_r_1_0_1/24 - H_r_1_1/12 + H_r_1_1_0/24 - H_r_1_1_1/24 - 31*pi**2/576 + 13/648 + zeta(3)/12
```

`eps^-1`:

```text
-7*pi**2*H_r_0/72 + H_r_0*zeta(3)/6 - 11*H_r_0_0/36 - pi**2*H_r_0_0/48 - 7*H_r_0_0_0/36 - H_r_0_0_0_0/12 + H_r_0_0_0_1/12 + 7*H_r_0_0_1/36 + H_r_0_0_1_0/12 - H_r_0_0_1_1/12 - pi**2*H_r_0_1/16 + 11*H_r_0_1/36 + 7*H_r_0_1_0/36 - H_r_0_1_0_0/12 + H_r_0_1_0_1/12 - 7*H_r_0_1_1/36 + H_r_0_1_1_0/12 - H_r_0_1_1_1/12 - 7*pi**2*H_r_1/72 + H_r_1*zeta(3)/6 - pi**2*H_r_1_0/16 + 11*H_r_1_0/36 - 7*H_r_1_0_0/36 - H_r_1_0_0_0/12 + H_r_1_0_0_1/12 + 7*H_r_1_0_1/36 + H_r_1_0_1_0/12 - H_r_1_0_1_1/12 - 11*H_r_1_1/36 - pi**2*H_r_1_1/48 + 7*H_r_1_1_0/36 - H_r_1_1_0_0/12 + H_r_1_1_0_1/12 - 7*H_r_1_1_1/36 + H_r_1_1_1_0/12 - H_r_1_1_1_1/12 - 245*pi**2/1728 - pi**4/240 + 101/486 + 49*zeta(3)/144
```

`eps^0`:

```text
-127*pi**2*H_r_0/432 - 19*pi**4*H_r_0/960 + 47*H_r_0*zeta(3)/72 - 47*H_r_0_0/54 - pi**2*H_r_0_0/18 + H_r_0_0*zeta(3)/12 - 127*H_r_0_0_0/216 + 7*pi**2*H_r_0_0_0/96 + 11*H_r_0_0_0_0/18 + 17*H_r_0_0_0_0_0/24 - 17*H_r_0_0_0_0_1/24 - 11*H_r_0_0_0_1/18 - 17*H_r_0_0_0_1_0/24 + 17*H_r_0_0_0_1_1/24 - 11*pi**2*H_r_0_0_1/96 + 127*H_r_0_0_1/216 - 11*H_r_0_0_1_0/18 - 19*H_r_0_0_1_0_0/24 + 19*H_r_0_0_1_0_1/24 + 11*H_r_0_0_1_1/18 + 19*H_r_0_0_1_1_0/24 - 19*H_r_0_0_1_1_1/24 - 13*pi**2*H_r_0_1/48 + 11*H_r_0_1*zeta(3)/24 + 47*H_r_0_1/54 - 11*pi**2*H_r_0_1_0/96 + 127*H_r_0_1_0/216 - 91*H_r_0_1_0_0/72 - 5*H_r_0_1_0_0_0/12 + 5*H_r_0_1_0_0_1/12 + 91*H_r_0_1_0_1/72 + 5*H_r_0_1_0_1_0/12 - 5*H_r_0_1_0_1_1/12 - 11*pi**2*H_r_0_1_1/96 - 127*H_r_0_1_1/216 + 91*H_r_0_1_1_0/72 - H_r_0_1_1_0_0/24 + H_r_0_1_1_0_1/24 - 91*H_r_0_1_1_1/72 + H_r_0_1_1_1_0/24 - H_r_0_1_1_1_1/24 - 127*pi**2*H_r_1/432 - 19*pi**4*H_r_1/960 + 47*H_r_1*zeta(3)/72 - 13*pi**2*H_r_1_0/48 + 11*H_r_1_0*zeta(3)/24 + 47*H_r_1_0/54 - 11*pi**2*H_r_1_0_0/96 - 127*H_r_1_0_0/216 - 91*H_r_1_0_0_0/72 - H_r_1_0_0_0_0/24 + H_r_1_0_0_0_1/24 + 91*H_r_1_0_0_1/72 + H_r_1_0_0_1_0/24 - H_r_1_0_0_1_1/24 - 11*pi**2*H_r_1_0_1/96 + 127*H_r_1_0_1/216 + 91*H_r_1_0_1_0/72 - 5*H_r_1_0_1_0_0/12 + 5*H_r_1_0_1_0_1/12 - 91*H_r_1_0_1_1/72 + 5*H_r_1_0_1_1_0/12 - 5*H_r_1_0_1_1_1/12 - 47*H_r_1_1/54 - pi**2*H_r_1_1/18 + H_r_1_1*zeta(3)/12 - 11*pi**2*H_r_1_1_0/96 + 127*H_r_1_1_0/216 + 11*H_r_1_1_0_0/18 - 19*H_r_1_1_0_0_0/24 + 19*H_r_1_1_0_0_1/24 - 11*H_r_1_1_0_1/18 + 19*H_r_1_1_0_1_0/24 - 19*H_r_1_1_0_1_1/24 - 127*H_r_1_1_1/216 + 7*pi**2*H_r_1_1_1/96 - 11*H_r_1_1_1_0/18 + 17*H_r_1_1_1_0_0/24 - 17*H_r_1_1_1_0_1/24 + 11*H_r_1_1_1_1/18 - 17*H_r_1_1_1_1_0/24 + 17*H_r_1_1_1_1_1/24 - 695*pi**4/13824 - 1217*pi**2/2592 + 182/729 + pi**2*zeta(3)/24 + zeta(5) + 599*zeta(3)/432
```

Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "expression_language": "SymPy; M_k_j means the original M(k,j,z), zeta_k=zeta(k)",
  "M_definition": "M(k,j,z)=integral_0^1 [log(beta*(1-beta))]^k [log(z)+log((1-beta)/beta)]^j/(beta-z/(1+z)) d beta, j>=1; removable internal point",
  "coefficient_field": "Q with explicit zeta constants and polynomial products of the original M at the same z; no additional rational z prefactor",
  "Hlog_expression_language": "Maple/HyperInt; zeta[n]=Zeta(n); Hlog uses dt/(t-letter)",
  "Hlog_variables": {
    "r": "1/(1+z)",
    "z": "positive real kinematic ratio; z1 in the raw source"
  },
  "Hlog_convention": {
    "empty_word": "1",
    "recursion": "Hlog(x,[a,w])=integral_0^x dt/(t-a) Hlog(t,w)",
    "all_zero_word": "log(x)^n/n!",
    "Hlog(x,[1])": "log(1-x)",
    "branch": "real tangential-basepoint GPL on 0<x<1"
  }
}
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-09/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-09/derivation.md).

## 3e-3g+Pi1-10

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-10/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-10/r003.json) (revision r003); SHA-256 `10d267285d6165845287d72f077c6c6384835348bafa7717367d0af869864413`.

Native unweighted field: `/laurent_Hlog_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — zero established by source pole metadata:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
1/768
```

`eps^-3` — explicit source coefficient:

```text
1/384 + pi**2/768
```

`eps^-2` — explicit source coefficient:

```text
-3*pi**2*H_z_m1/256 - 3*H_z_m1_0_0/128 + rho**3/512 + 3*pi**2*rho/512 - 7*zeta(3)/384 + 1/192 + 7*pi**2/3072
```

`eps^-1` — explicit source coefficient:

```text
-19*pi**2*H_z_0_m1/512 - 19*H_z_0_m1_0_0/256 - H_z_m1*(-19*zeta(3)/128 + 3*pi**2/128) - 19*pi**2*H_z_m1_0/512 - 3*H_z_m1_0_0/64 - 19*H_z_m1_0_0_0/256 + 19*pi**2*H_z_m1_m1/256 + 19*H_z_m1_m1_0_0/128 + 13*rho**4/6144 + rho**3/256 + 23*pi**2*rho**2/3072 - rho*(-3*pi**2/256 + 19*zeta(3)/256) - 11*zeta(3)/256 + 7/384 + 11*pi**2/1536 + 43*pi**4/6144
```

`eps^0` — explicit source coefficient:

```text
-17*pi**2*H_z_0_0_m1/256 - 17*H_z_0_0_m1_0_0/128 - H_z_0_m1*(-81*zeta(3)/256 + 19*pi**2/256) - 11*pi**2*H_z_0_m1_0/192 - 19*H_z_0_m1_0_0/128 - 7*H_z_0_m1_0_0_0/256 + 81*pi**2*H_z_0_m1_m1/512 + 81*H_z_0_m1_m1_0_0/256 - H_z_m1*(-19*zeta(3)/64 + 3*pi**2/64 + 683*pi**4/30720) - H_z_m1_0*(-81*zeta(3)/256 + 19*pi**2/256) - H_z_m1_0_0*(3/32 + 41*pi**2/768) - 19*H_z_m1_0_0_0/128 + 5*H_z_m1_0_0_0_0/256 + 81*pi**2*H_z_m1_0_m1/512 + 81*H_z_m1_0_m1_0_0/256 - H_z_m1_m1*(-19*pi**2/128 + 81*zeta(3)/128) + 81*pi**2*H_z_m1_m1_0/512 + 19*H_z_m1_m1_0_0/64 + 81*H_z_m1_m1_0_0_0/256 - 81*pi**2*H_z_m1_m1_m1/256 - 81*H_z_m1_m1_m1_0_0/128 - 3*rho**5/4096 + 13*rho**4/3072 - rho**3*(-29*pi**2/3072 - 3/64)/6 - rho**2*(-23*pi**2/768 + 47*zeta(3)/256)/2 - rho*(-2287*pi**4/184320 - 3*pi**2/128 + 19*zeta(3)/128) - 49*pi**2*zeta(3)/384 - 521*zeta(5)/768 - 7*zeta(3)/128 + 7/192 + pi**2/128 + 1657*pi**4/122880
```

### Af contribution

Consumed field: `/laurent_Hlog_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
-1/288
```

`eps^-3`:

```text
-pi**2/288 - 5/864
```

`eps^-2`:

```text
pi**2*H_z_m1/32 + H_z_m1_0_0/16 - rho**3/192 - pi**2*rho/64 - 17*pi**2/3456 - 7/648 + 7*zeta(3)/144
```

`eps^-1`:

```text
19*pi**2*H_z_0_m1/192 + 19*H_z_0_m1_0_0/96 - 19*H_z_m1*zeta(3)/48 + 5*pi**2*H_z_m1/96 + 19*pi**2*H_z_m1_0/192 + 5*H_z_m1_0_0/48 + 19*H_z_m1_0_0_0/96 - 19*pi**2*H_z_m1_m1/96 - 19*H_z_m1_m1_0_0/48 - 13*rho**4/2304 - 5*rho**3/576 - 23*pi**2*rho**2/1152 - 5*pi**2*rho/192 + 19*rho*zeta(3)/96 - 43*pi**4/2304 - 169*pi**2/10368 - 163/3888 + 85*zeta(3)/864
```

`eps^0`:

```text
17*pi**2*H_z_0_0_m1/96 + 17*H_z_0_0_m1_0_0/48 - 27*H_z_0_m1*zeta(3)/32 + 95*pi**2*H_z_0_m1/576 + 11*pi**2*H_z_0_m1_0/72 + 95*H_z_0_m1_0_0/288 + 7*H_z_0_m1_0_0_0/96 - 27*pi**2*H_z_0_m1_m1/64 - 27*H_z_0_m1_m1_0_0/32 - 95*H_z_m1*zeta(3)/144 + 7*pi**2*H_z_m1/72 + 683*pi**4*H_z_m1/11520 - 27*H_z_m1_0*zeta(3)/32 + 95*pi**2*H_z_m1_0/576 + 7*H_z_m1_0_0/36 + 41*pi**2*H_z_m1_0_0/288 + 95*H_z_m1_0_0_0/288 - 5*H_z_m1_0_0_0_0/96 - 27*pi**2*H_z_m1_0_m1/64 - 27*H_z_m1_0_m1_0_0/32 - 95*pi**2*H_z_m1_m1/288 + 27*H_z_m1_m1*zeta(3)/16 - 27*pi**2*H_z_m1_m1_0/64 - 95*H_z_m1_m1_0_0/144 - 27*H_z_m1_m1_0_0_0/32 + 27*pi**2*H_z_m1_m1_m1/32 + 27*H_z_m1_m1_m1_0_0/16 + rho**5/512 - 65*rho**4/6912 - 29*pi**2*rho**3/6912 - 7*rho**3/432 - 115*pi**2*rho**2/3456 + 47*rho**2*zeta(3)/192 - 2287*pi**4*rho/69120 - 7*pi**2*rho/144 + 95*rho*zeta(3)/288 - 4111*pi**4/138240 - 49*pi**2/3888 - 893/11664 + 251*zeta(3)/2592 + 521*zeta(5)/288 + 49*pi**2*zeta(3)/144
```

Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "expression_language": "SymPy; original M_k_j, N_k_j at z; zeta_k=zeta(k)",
  "M_definition": "M(k,j,z)=integral_0^1 [log(beta*(1-beta))]^k [log(z)+log((1-beta)/beta)]^j/(beta-z/(1+z)) d beta, j>=1; removable internal point",
  "N_definition": "N(k,j,z)=M(k,j+1,z)-log(z)*M(k,j,z)",
  "coefficient_field": "Q with explicit zeta constants and polynomial products of the original M/N at the same z",
  "Hlog_language": "Maple/HyperInt; Hlog(z,word) with letters 0,-1; Hlog(z,[0])=log(z), Hlog(z,[-1])=log(1+z)"
}
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-10/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-10/derivation.md).

## 3e-3g+Pi1-11

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-11/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-11/r003.json) (revision r003); SHA-256 `2617946d55e88e9dcbcf1457605d332fc2ee5170165a930399ad1489c3f00c35`.

Native unweighted field: `/laurent_Hlog_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — zero established by source pole metadata:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
-1/384 + pi**2/1152
```

`eps^-2` — explicit source coefficient:

```text
-pi**2*H_z_m1/128 - H_z_m1_0_0/64 + rho**3/768 - 3*rho**2/512 + pi**2*rho/256 - 19*pi**2/4608 - zeta(3)/96 - 1/384
```

`eps^-1` — explicit source coefficient:

```text
-7*pi**2*H_z_0_m1/256 - 7*H_z_0_m1_0_0/128 - H_z_m1*(-7*zeta(3)/64 - pi**2/128) - 7*pi**2*H_z_m1_0/256 + H_z_m1_0_0/64 - 7*H_z_m1_0_0_0/128 + 7*pi**2*H_z_m1_m1/128 + 7*H_z_m1_m1_0_0/64 + rho**4/512 - rho**3*(3*z/512 + 5/512)/(6*(z + 1)) - rho**2*(3/256 - 5*pi**2/384)/2 - rho*(11*pi**2*z/3072 + 7*z*zeta(3)/128 + 13*pi**2/3072 + 7*zeta(3)/128)/(z + 1) - pi**2/144 - 3/128 + 5*zeta(3)/192 + 89*pi**4/13824
```

`eps^0` — explicit source coefficient:

```text
-7*pi**2*H_z_0_0_m1/128 - 7*H_z_0_0_m1_0_0/64 - H_z_0_m1*(-27*zeta(3)/128 + 5*pi**2/256) - pi**2*H_z_0_m1_0/24 - 5*H_z_0_m1_0_0/128 - 5*H_z_0_m1_0_0_0/128 + 27*pi**2*H_z_0_m1_m1/256 + 27*H_z_0_m1_m1_0_0/128 - H_z_m1*(-5*zeta(3)/64 + pi**2/128 + 103*pi**4/7680) - H_z_m1_0*(-27*zeta(3)/128 + 5*pi**2/256) - H_z_m1_0_0*(1/64 + 25*pi**2/768) - 5*H_z_m1_0_0_0/128 + H_z_m1_0_0_0_0/64 + 27*pi**2*H_z_m1_0_m1/256 + 27*H_z_m1_0_m1_0_0/128 - H_z_m1_m1*(-5*pi**2/128 + 27*zeta(3)/64) + 27*pi**2*H_z_m1_m1_0/256 + 5*H_z_m1_m1_0_0/64 + 27*H_z_m1_m1_0_0_0/128 - 27*pi**2*H_z_m1_m1_m1/128 - 27*H_z_m1_m1_m1_0_0/64 - rho**5/3072 + 15*rho**4/2048 + rho**3*(5*z/512 + 3*pi**2*z/256 + 3/512 + 3*pi**2/256)/(6*(z + 1)) - rho**2*(-59*pi**2/3072 + 9/256 + 11*zeta(3)/64)/2 + rho*(-5*z*zeta(3)/128 + 13*pi**2*z/3072 + 197*pi**4*z/23040 - 5*zeta(3)/128 + 11*pi**2/3072 + 197*pi**4/23040)/(z + 1) - 41*pi**2*zeta(3)/288 - 5*zeta(5)/12 - 17*pi**2/2304 - 17*zeta(3)/384 - 17/384 + 997*pi**4/138240
```

### Af contribution

Consumed field: `/laurent_Hlog_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
1/144 - pi**2/432
```

`eps^-2`:

```text
pi**2*H_z_m1/48 + H_z_m1_0_0/24 - rho**3/288 + rho**2/64 - pi**2*rho/96 + 1/216 + zeta(3)/36 + 61*pi**2/5184
```

`eps^-1`:

```text
(2268*pi**2*H_z_0_m1*z + 2268*pi**2*H_z_0_m1 + 4536*H_z_0_m1_0_0*z + 4536*H_z_0_m1_0_0 - 9072*H_z_m1*z*zeta(3) - 864*pi**2*H_z_m1*z - 9072*H_z_m1*zeta(3) - 864*pi**2*H_z_m1 + 2268*pi**2*H_z_m1_0*z + 2268*pi**2*H_z_m1_0 - 1728*H_z_m1_0_0*z - 1728*H_z_m1_0_0 + 4536*H_z_m1_0_0_0*z + 4536*H_z_m1_0_0_0 - 4536*pi**2*H_z_m1_m1*z - 4536*pi**2*H_z_m1_m1 - 9072*H_z_m1_m1_0_0*z - 9072*H_z_m1_m1_0_0 - 162*rho**4*z - 162*rho**4 + 117*rho**3*z + 171*rho**3 - 540*pi**2*rho**2*z + 324*rho**2*z - 540*pi**2*rho**2 + 324*rho**2 + 405*pi**2*rho*z + 4536*rho*z*zeta(3) + 459*pi**2*rho + 4536*rho*zeta(3) - 534*pi**4*z - 2448*z*zeta(3) + 1824*z + 478*pi**2*z - 534*pi**4 - 2448*zeta(3) + 1824 + 478*pi**2)/(31104*z + 31104)
```

`eps^0`:

```text
(136080*pi**2*H_z_0_0_m1*z + 136080*pi**2*H_z_0_0_m1 + 272160*H_z_0_0_m1_0_0*z + 272160*H_z_0_0_m1_0_0 - 524880*H_z_0_m1*z*zeta(3) + 25920*pi**2*H_z_0_m1*z - 524880*H_z_0_m1*zeta(3) + 25920*pi**2*H_z_0_m1 + 103680*pi**2*H_z_0_m1_0*z + 103680*pi**2*H_z_0_m1_0 + 51840*H_z_0_m1_0_0*z + 51840*H_z_0_m1_0_0 + 97200*H_z_0_m1_0_0_0*z + 97200*H_z_0_m1_0_0_0 - 262440*pi**2*H_z_0_m1_m1*z - 262440*pi**2*H_z_0_m1_m1 - 524880*H_z_0_m1_m1_0_0*z - 524880*H_z_0_m1_m1_0_0 - 103680*H_z_m1*z*zeta(3) + 21600*pi**2*H_z_m1*z + 33372*pi**4*H_z_m1*z - 103680*H_z_m1*zeta(3) + 21600*pi**2*H_z_m1 + 33372*pi**4*H_z_m1 - 524880*H_z_m1_0*z*zeta(3) + 25920*pi**2*H_z_m1_0*z - 524880*H_z_m1_0*zeta(3) + 25920*pi**2*H_z_m1_0 + 43200*H_z_m1_0_0*z + 81000*pi**2*H_z_m1_0_0*z + 43200*H_z_m1_0_0 + 81000*pi**2*H_z_m1_0_0 + 51840*H_z_m1_0_0_0*z + 51840*H_z_m1_0_0_0 - 38880*H_z_m1_0_0_0_0*z - 38880*H_z_m1_0_0_0_0 - 262440*pi**2*H_z_m1_0_m1*z - 262440*pi**2*H_z_m1_0_m1 - 524880*H_z_m1_0_m1_0_0*z - 524880*H_z_m1_0_m1_0_0 - 51840*pi**2*H_z_m1_m1*z + 1049760*H_z_m1_m1*z*zeta(3) - 51840*pi**2*H_z_m1_m1 + 1049760*H_z_m1_m1*zeta(3) - 262440*pi**2*H_z_m1_m1_0*z - 262440*pi**2*H_z_m1_m1_0 - 103680*H_z_m1_m1_0_0*z - 103680*H_z_m1_m1_0_0 - 524880*H_z_m1_m1_0_0_0*z - 524880*H_z_m1_m1_0_0_0 + 524880*pi**2*H_z_m1_m1_m1*z + 524880*pi**2*H_z_m1_m1_m1 + 1049760*H_z_m1_m1_m1_0_0*z + 1049760*H_z_m1_m1_m1_0_0 + 810*rho**5*z + 810*rho**5 - 16605*rho**4*z - 16605*rho**4 - 4860*pi**2*rho**3*z - 4140*rho**3*z - 4860*pi**2*rho**3 - 3060*rho**3 - 18495*pi**2*rho**2*z + 35640*rho**2*z + 213840*rho**2*z*zeta(3) - 18495*pi**2*rho**2 + 35640*rho**2 + 213840*rho**2*zeta(3) - 21276*pi**4*rho*z - 11340*pi**2*rho*z + 51840*rho*z*zeta(3) - 21276*pi**4*rho - 10260*pi**2*rho + 51840*rho*zeta(3) - 12606*pi**4*z + 88320*z + 10640*pi**2*z + 126000*z*zeta(3) + 1036800*z*zeta(5) + 354240*pi**2*z*zeta(3) - 12606*pi**4 + 88320 + 10640*pi**2 + 126000*zeta(3) + 1036800*zeta(5) + 354240*pi**2*zeta(3))/(933120*z + 933120)
```

Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "expression_language": "SymPy; original M_k_j, N_k_j at z; zeta_k=zeta(k); chi=(1-z)/(1+z)",
  "M_definition": "M(k,j,z)=integral_0^1 [log(beta*(1-beta))]^k [log(z)+log((1-beta)/beta)]^j/(beta-z/(1+z)) d beta, j>=1; removable internal point",
  "N_definition": "N(k,j,z)=M(k,j+1,z)-log(z)*M(k,j,z)",
  "chi_definition": "chi=(1-z)/(1+z)",
  "coefficient_field": "Q(z) with explicit zeta constants and polynomial products of the original M/N at the same z",
  "Hlog_language": "Maple/HyperInt; Hlog(z,word) with letters 0,-1; Hlog(z,[0])=log(z), Hlog(z,[-1])=log(1+z)"
}
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-11/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-11/derivation.md).

## 3e-3g+Pi1-12

Source: [benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-12/r003.json](sources/benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-12/r003.json) (revision r003); SHA-256 `f81e2f392727a44c8016606af30e394080e9668f35be94373332054b8976a92b`.

Native unweighted field: `/laurent_Hlog_coefficients`. native scalar I before the external matter coefficient Bf(eps).

### Native unweighted coefficients

`eps^-6` — zero established by source pole metadata:

```text
0
```

`eps^-5` — zero established by source pole metadata:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
1/768
```

`eps^-3` — explicit source coefficient:

```text
1/384 + pi**2/768
```

`eps^-2` — explicit source coefficient:

```text
-3*pi**2*H_z_m1/256 - 3*H_z_m1_0_0/128 + rho**3/512 + 3*pi**2*rho/512 - 7*zeta(3)/384 + 1/192 + 7*pi**2/3072
```

`eps^-1` — explicit source coefficient:

```text
-19*pi**2*H_z_0_m1/512 - 19*H_z_0_m1_0_0/256 - H_z_m1*(-19*zeta(3)/128 + 3*pi**2/128) - 19*pi**2*H_z_m1_0/512 - 3*H_z_m1_0_0/64 - 19*H_z_m1_0_0_0/256 + 19*pi**2*H_z_m1_m1/256 + 19*H_z_m1_m1_0_0/128 + 13*rho**4/6144 + rho**3/256 + 23*pi**2*rho**2/3072 - rho*(-3*pi**2/256 + 19*zeta(3)/256) - 11*zeta(3)/256 + 7/384 + 11*pi**2/1536 + 43*pi**4/6144
```

`eps^0` — explicit source coefficient:

```text
-17*pi**2*H_z_0_0_m1/256 - 17*H_z_0_0_m1_0_0/128 - H_z_0_m1*(-81*zeta(3)/256 + 19*pi**2/256) - 11*pi**2*H_z_0_m1_0/192 - 19*H_z_0_m1_0_0/128 - 7*H_z_0_m1_0_0_0/256 + 81*pi**2*H_z_0_m1_m1/512 + 81*H_z_0_m1_m1_0_0/256 - H_z_m1*(-19*zeta(3)/64 + 3*pi**2/64 + 683*pi**4/30720) - H_z_m1_0*(-81*zeta(3)/256 + 19*pi**2/256) - H_z_m1_0_0*(3/32 + 41*pi**2/768) - 19*H_z_m1_0_0_0/128 + 5*H_z_m1_0_0_0_0/256 + 81*pi**2*H_z_m1_0_m1/512 + 81*H_z_m1_0_m1_0_0/256 - H_z_m1_m1*(-19*pi**2/128 + 81*zeta(3)/128) + 81*pi**2*H_z_m1_m1_0/512 + 19*H_z_m1_m1_0_0/64 + 81*H_z_m1_m1_0_0_0/256 - 81*pi**2*H_z_m1_m1_m1/256 - 81*H_z_m1_m1_m1_0_0/128 - 3*rho**5/4096 + 13*rho**4/3072 - rho**3*(-29*pi**2/3072 - 3/64)/6 - rho**2*(-23*pi**2/768 + 47*zeta(3)/256)/2 - rho*(-2287*pi**4/184320 - 3*pi**2/128 + 19*zeta(3)/128) - 49*pi**2*zeta(3)/384 - 521*zeta(5)/768 - 7*zeta(3)/128 + 7/192 + pi**2/128 + 1657*pi**4/122880
```

### Af contribution

Consumed field: `/laurent_Hlog_coefficients`. External scalar multiplier: `(8*eps - 8)/(3 - 2*eps)`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
-1/288
```

`eps^-3`:

```text
-pi**2/288 - 5/864
```

`eps^-2`:

```text
pi**2*H_z_m1/32 + H_z_m1_0_0/16 - rho**3/192 - pi**2*rho/64 - 17*pi**2/3456 - 7/648 + 7*zeta(3)/144
```

`eps^-1`:

```text
19*pi**2*H_z_0_m1/192 + 19*H_z_0_m1_0_0/96 - 19*H_z_m1*zeta(3)/48 + 5*pi**2*H_z_m1/96 + 19*pi**2*H_z_m1_0/192 + 5*H_z_m1_0_0/48 + 19*H_z_m1_0_0_0/96 - 19*pi**2*H_z_m1_m1/96 - 19*H_z_m1_m1_0_0/48 - 13*rho**4/2304 - 5*rho**3/576 - 23*pi**2*rho**2/1152 - 5*pi**2*rho/192 + 19*rho*zeta(3)/96 - 43*pi**4/2304 - 169*pi**2/10368 - 163/3888 + 85*zeta(3)/864
```

`eps^0`:

```text
17*pi**2*H_z_0_0_m1/96 + 17*H_z_0_0_m1_0_0/48 - 27*H_z_0_m1*zeta(3)/32 + 95*pi**2*H_z_0_m1/576 + 11*pi**2*H_z_0_m1_0/72 + 95*H_z_0_m1_0_0/288 + 7*H_z_0_m1_0_0_0/96 - 27*pi**2*H_z_0_m1_m1/64 - 27*H_z_0_m1_m1_0_0/32 - 95*H_z_m1*zeta(3)/144 + 7*pi**2*H_z_m1/72 + 683*pi**4*H_z_m1/11520 - 27*H_z_m1_0*zeta(3)/32 + 95*pi**2*H_z_m1_0/576 + 7*H_z_m1_0_0/36 + 41*pi**2*H_z_m1_0_0/288 + 95*H_z_m1_0_0_0/288 - 5*H_z_m1_0_0_0_0/96 - 27*pi**2*H_z_m1_0_m1/64 - 27*H_z_m1_0_m1_0_0/32 - 95*pi**2*H_z_m1_m1/288 + 27*H_z_m1_m1*zeta(3)/16 - 27*pi**2*H_z_m1_m1_0/64 - 95*H_z_m1_m1_0_0/144 - 27*H_z_m1_m1_0_0_0/32 + 27*pi**2*H_z_m1_m1_m1/32 + 27*H_z_m1_m1_m1_0_0/16 + rho**5/512 - 65*rho**4/6912 - 29*pi**2*rho**3/6912 - 7*rho**3/432 - 115*pi**2*rho**2/3456 + 47*rho**2*zeta(3)/192 - 2287*pi**4*rho/69120 - 7*pi**2*rho/144 + 95*rho*zeta(3)/288 - 4111*pi**4/138240 - 49*pi**2/3888 - 893/11664 + 251*zeta(3)/2592 + 521*zeta(5)/288 + 49*pi**2*zeta(3)/144
```

Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "expression_language": "SymPy; original M_k_j, N_k_j at z; zeta_k=zeta(k)",
  "M_definition": "M(k,j,z)=integral_0^1 [log(beta*(1-beta))]^k [log(z)+log((1-beta)/beta)]^j/(beta-z/(1+z)) d beta, j>=1; removable internal point",
  "N_definition": "N(k,j,z)=M(k,j+1,z)-log(z)*M(k,j,z)",
  "coefficient_field": "Q with explicit zeta constants and polynomial products of the original M/N at the same z",
  "Hlog_language": "Maple/HyperInt; Hlog(z,word) with letters 0,-1; Hlog(z,[0])=log(z), Hlog(z,[-1])=log(1+z)"
}
```

Source dependency: [workflow/checks/tgv-pi1-representation-r003/verification.json](sources/workflow/checks/tgv-pi1-representation-r003/verification.json).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-12/derivation.md](sources/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-12/derivation.md).

## 3e-Gamma1-01

Source: [benchmarks/self-energy/gamma1/results/3e-Gamma1-01/r002.json](sources/benchmarks/self-energy/gamma1/results/3e-Gamma1-01/r002.json) (revision r002); SHA-256 `414c5ca5942257846ae88b083d9b4c8c6a135370dbff8f49563a104568c20c33`.

Native unweighted field: `/channels/q/I_coefficients`. quark-channel scalar I_q after stripping the external color monomial.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/864
```

`eps^-4` — explicit source coefficient:

```text
-7/2592
```

`eps^-3` — explicit source coefficient:

```text
-131/7776 + 19*pi**2/10368
```

`eps^-2` — explicit source coefficient:

```text
-3251/46656 + 59*pi**2/31104 + 17*zeta(3)/288
```

`eps^-1` — explicit source coefficient:

```text
-68435/279936 - 47*zeta(3)/864 + 127*pi**2/93312 + 827*pi**4/414720
```

`eps^0` — explicit source coefficient:

```text
-1365323/1679616 - 305*zeta(3)/2592 - 2921*pi**2/559872 + 739*pi**4/1244160 + 11*pi**2*zeta(3)/384 + 1331*zeta(5)/1440
```

### Af contribution

Consumed field: `/channels/q/I_coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
1/864
```

`eps^-4`:

```text
-7/2592
```

`eps^-3`:

```text
-131/7776 + 19*pi**2/10368
```

`eps^-2`:

```text
-3251/46656 + 59*pi**2/31104 + 17*zeta(3)/288
```

`eps^-1`:

```text
-68435/279936 - 47*zeta(3)/864 + 127*pi**2/93312 + 827*pi**4/414720
```

`eps^0`:

```text
-1365323/1679616 - 305*zeta(3)/2592 - 2921*pi**2/559872 + 739*pi**4/1244160 + 11*pi**2*zeta(3)/384 + 1331*zeta(5)/1440
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "coefficient_notation": {
    "Z2": "zeta(2)",
    "Z3": "zeta(3)",
    "Z5": "zeta(5)"
  }
}
```

Source dependency: [benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-01/derivation.md](sources/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-01/derivation.md).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-01/derivation.md](sources/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-01/derivation.md).

## 3e-Gamma1-02

Source: [benchmarks/self-energy/gamma1/results/3e-Gamma1-02/r002.json](sources/benchmarks/self-energy/gamma1/results/3e-Gamma1-02/r002.json) (revision r002); SHA-256 `7c77526e84b52809aab0990fde75cf37789978741c0af92696ffb339c31f6a8b`.

Native unweighted field: `/channels/q/I_coefficients`. quark-channel scalar I_q after stripping the external color monomial.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
1/864
```

`eps^-4` — explicit source coefficient:

```text
-7/2592
```

`eps^-3` — explicit source coefficient:

```text
-131/7776 + 19*pi**2/10368
```

`eps^-2` — explicit source coefficient:

```text
-3251/46656 + 59*pi**2/31104 + 17*zeta(3)/288
```

`eps^-1` — explicit source coefficient:

```text
-68435/279936 - 47*zeta(3)/864 + 127*pi**2/93312 + 827*pi**4/414720
```

`eps^0` — explicit source coefficient:

```text
-1365323/1679616 - 305*zeta(3)/2592 - 2921*pi**2/559872 + 739*pi**4/1244160 + 11*pi**2*zeta(3)/384 + 1331*zeta(5)/1440
```

### Af contribution

Consumed field: `/channels/q/I_coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
1/864
```

`eps^-4`:

```text
-7/2592
```

`eps^-3`:

```text
-131/7776 + 19*pi**2/10368
```

`eps^-2`:

```text
-3251/46656 + 59*pi**2/31104 + 17*zeta(3)/288
```

`eps^-1`:

```text
-68435/279936 - 47*zeta(3)/864 + 127*pi**2/93312 + 827*pi**4/414720
```

`eps^0`:

```text
-1365323/1679616 - 305*zeta(3)/2592 - 2921*pi**2/559872 + 739*pi**4/1244160 + 11*pi**2*zeta(3)/384 + 1331*zeta(5)/1440
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "coefficient_notation": {
    "Z2": "zeta(2)",
    "Z3": "zeta(3)",
    "Z5": "zeta(5)"
  }
}
```

Source dependency: [benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-02/derivation.md](sources/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-02/derivation.md).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-02/derivation.md](sources/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-02/derivation.md).

## 3e-Gamma1-03

Source: [benchmarks/self-energy/gamma1/results/3e-Gamma1-03/z-r003.json](sources/benchmarks/self-energy/gamma1/results/3e-Gamma1-03/z-r003.json) (revision z-r003); SHA-256 `68afb9c88590c00b58a8dea4fba66f4f820416be9c139d5a5dda556f8daa4e17`.

Native unweighted field: `/channels/q/I_coefficients`. quark-channel scalar I_q after stripping the external color monomial.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
0
```

`eps^-2` — explicit source coefficient:

```text
-pi**2*H_z_m1/48 - H_z_m1_0_0/24 + rho**3/288 + rho**2/48 + pi**2*rho/96 - zeta(3)/24 + pi**2/54
```

`eps^-1` — explicit source coefficient:

```text
-pi**2*H_z_0_m1/32 - H_z_0_m1_0_0/16 - H_z_m1*(-zeta(3)/8 + 25*pi**2/288) - pi**2*H_z_m1_0/32 - 25*H_z_m1_0_0/144 - H_z_m1_0_0_0/16 + pi**2*H_z_m1_m1/16 + H_z_m1_m1_0_0/8 + rho**4/576 + rho**3*(11*z/144 + 7/72)/(6*(z + 1)) - rho**2*(-53/288 - pi**2/144)/2 + rho*(-z*zeta(3)/16 + pi**2*z/24 - zeta(3)/16 + 13*pi**2/288)/(z + 1) - zeta(3)/16 - 1/72 + pi**4/1440 + 445*pi**2/5184
```

`eps^0` — explicit source coefficient:

```text
pi**2*H_z_0_0_m1/24 + H_z_0_0_m1_0_0/12 - H_z_0_m1*(-3*zeta(3)/16 + 25*pi**2/192) - pi**2*H_z_0_m1_0/72 - 25*H_z_0_m1_0_0/96 + 5*H_z_0_m1_0_0_0/48 + 3*pi**2*H_z_0_m1_m1/32 + 3*H_z_0_m1_m1_0_0/16 - H_z_m1*(-25*zeta(3)/48 + 7*pi**4/360 + 577*pi**2/1728) - H_z_m1_0*(-3*zeta(3)/16 + 25*pi**2/192) - H_z_m1_0_0*(5*pi**2/288 + 577/864) - 25*H_z_m1_0_0_0/96 + H_z_m1_0_0_0_0/3 + 3*pi**2*H_z_m1_0_m1/32 + 3*H_z_m1_0_m1_0_0/16 - H_z_m1_m1*(-25*pi**2/96 + 3*zeta(3)/8) + 3*pi**2*H_z_m1_m1_0/32 + 25*H_z_m1_m1_0_0/48 + 3*H_z_m1_m1_0_0_0/16 - 3*pi**2*H_z_m1_m1_m1/16 - 3*H_z_m1_m1_m1_0_0/8 - 17*rho**5/5760 - 11*rho**4/864 - rho**3*(-565*z/1728 + 5*pi**2*z/96 - 589/1728 + 5*pi**2/96)/(6*(z + 1)) - rho**2*(-1193/1728 - 19*pi**2/864 + 3*zeta(3)/8)/2 + rho*(-25*z*zeta(3)/96 + pi**4*z/540 + 191*pi**2*z/1152 - 25*zeta(3)/96 + pi**4/540 + 581*pi**2/3456)/(z + 1) - 23*pi**2*zeta(3)/96 - 283*zeta(3)/288 - 1/16 + 5*zeta(5)/12 + 323*pi**4/17280 + 10609*pi**2/31104
```

### Af contribution

Consumed field: `/channels/q/I_coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
0
```

`eps^-2`:

```text
-pi**2*H_z_m1/48 - H_z_m1_0_0/24 + rho**3/288 + rho**2/48 + pi**2*rho/96 - zeta(3)/24 + pi**2/54
```

`eps^-1`:

```text
-pi**2*H_z_0_m1/32 - H_z_0_m1_0_0/16 - H_z_m1*(-zeta(3)/8 + 25*pi**2/288) - pi**2*H_z_m1_0/32 - 25*H_z_m1_0_0/144 - H_z_m1_0_0_0/16 + pi**2*H_z_m1_m1/16 + H_z_m1_m1_0_0/8 + rho**4/576 + rho**3*(11*z/144 + 7/72)/(6*z + 6) - rho**2*(-53/288 - pi**2/144)/2 + rho*(-z*zeta(3)/16 + pi**2*z/24 - zeta(3)/16 + 13*pi**2/288)/(z + 1) - zeta(3)/16 - 1/72 + pi**4/1440 + 445*pi**2/5184
```

`eps^0`:

```text
pi**2*H_z_0_0_m1/24 + H_z_0_0_m1_0_0/12 - H_z_0_m1*(-3*zeta(3)/16 + 25*pi**2/192) - pi**2*H_z_0_m1_0/72 - 25*H_z_0_m1_0_0/96 + 5*H_z_0_m1_0_0_0/48 + 3*pi**2*H_z_0_m1_m1/32 + 3*H_z_0_m1_m1_0_0/16 - H_z_m1*(-25*zeta(3)/48 + 7*pi**4/360 + 577*pi**2/1728) - H_z_m1_0*(-3*zeta(3)/16 + 25*pi**2/192) - H_z_m1_0_0*(5*pi**2/288 + 577/864) - 25*H_z_m1_0_0_0/96 + H_z_m1_0_0_0_0/3 + 3*pi**2*H_z_m1_0_m1/32 + 3*H_z_m1_0_m1_0_0/16 - H_z_m1_m1*(-25*pi**2/96 + 3*zeta(3)/8) + 3*pi**2*H_z_m1_m1_0/32 + 25*H_z_m1_m1_0_0/48 + 3*H_z_m1_m1_0_0_0/16 - 3*pi**2*H_z_m1_m1_m1/16 - 3*H_z_m1_m1_m1_0_0/8 - 17*rho**5/5760 - 11*rho**4/864 - rho**3*(-565*z/1728 + 5*pi**2*z/96 - 589/1728 + 5*pi**2/96)/(6*z + 6) - rho**2*(-1193/1728 - 19*pi**2/864 + 3*zeta(3)/8)/2 + rho*(-25*z*zeta(3)/96 + pi**4*z/540 + 191*pi**2*z/1152 - 25*zeta(3)/96 + pi**4/540 + 581*pi**2/3456)/(z + 1) - 23*pi**2*zeta(3)/96 - 283*zeta(3)/288 - 1/16 + 5*zeta(5)/12 + 323*pi**4/17280 + 10609*pi**2/31104
```

Existing compact or alternate representations are attached in the JSON with their source definitions and recorded validation. They were not converted anew for this export.

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "coefficient_notation": {
    "language": "Maple exact rational expressions",
    "argument": "z",
    "alphabet": [
      0,
      -1
    ],
    "Hlog(z,[a1,...,an])": "G(a1,...,an;z), kernels dt/(t-a1)",
    "Hlog(z,[-1])": "log(1+z)",
    "Hlog(z,[0])": "log(z)",
    "zeta[k]": "Riemann zeta(k)",
    "^": "power"
  }
}
```

Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-03/derivation.md](sources/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-03/derivation.md).

## 3e-Gamma1-04

Source: [benchmarks/self-energy/gamma1/results/3e-Gamma1-04/r003.json](sources/benchmarks/self-energy/gamma1/results/3e-Gamma1-04/r003.json) (revision r003); SHA-256 `473a7c442c48c1da04a2fb9287b6b02becd2c254f6cb7045a653c57a5c2d12c8`.

Native unweighted field: `/channels/q/I_coefficients`. quark-channel scalar I_q after stripping the external color monomial.

### Native unweighted coefficients

`eps^-6` — explicit source coefficient:

```text
0
```

`eps^-5` — explicit source coefficient:

```text
0
```

`eps^-4` — explicit source coefficient:

```text
0
```

`eps^-3` — explicit source coefficient:

```text
pi**2/432
```

`eps^-2` — explicit source coefficient:

```text
-pi**2*H_z_m1/48 - H_z_m1_0_0/24 + rho**3/288 + pi**2*rho/96 - 5*zeta(3)/72 + pi**2/162
```

`eps^-1` — explicit source coefficient:

```text
-3*pi**2*H_z_0_m1/32 - 3*H_z_0_m1_0_0/16 - H_z_m1*(-3*zeta(3)/8 + 13*pi**2/288) - 3*pi**2*H_z_m1_0/32 - 13*H_z_m1_0_0/144 - 3*H_z_m1_0_0_0/16 + 3*pi**2*H_z_m1_m1/16 + 3*H_z_m1_m1_0_0/8 + rho**4/192 - rho**3*(-z/18 - 5/144)/(6*(z + 1)) + pi**2*rho**2/48 + rho*(-3*z*zeta(3)/16 + 7*pi**2*z/288 - 3*zeta(3)/16 + pi**2/48)/(z + 1) - 83*zeta(3)/432 + 13*pi**2/972 + 11*pi**4/540
```

`eps^0` — explicit source coefficient:

```text
-13*pi**2*H_z_0_0_m1/48 - 13*H_z_0_0_m1_0_0/24 - H_z_0_m1*(-41*zeta(3)/48 + 11*pi**2/64) - 25*pi**2*H_z_0_m1_0/144 - 11*H_z_0_m1_0_0/32 - 3*H_z_0_m1_0_0_0/16 + 41*pi**2*H_z_0_m1_m1/96 + 41*H_z_0_m1_m1_0_0/48 - H_z_m1*(-11*zeta(3)/16 + 31*pi**2/432 + 77*pi**4/1440) - H_z_m1_0*(-41*zeta(3)/48 + 11*pi**2/64) - H_z_m1_0_0*(31/216 + 43*pi**2/288) - 11*H_z_m1_0_0_0/32 - H_z_m1_0_0_0_0/6 + 41*pi**2*H_z_m1_0_m1/96 + 41*H_z_m1_0_m1_0_0/48 - H_z_m1_m1*(-11*pi**2/32 + 41*zeta(3)/24) + 41*pi**2*H_z_m1_m1_0/96 + 11*H_z_m1_m1_0_0/16 + 41*H_z_m1_m1_0_0_0/48 - 41*pi**2*H_z_m1_m1_m1/48 - 41*H_z_m1_m1_m1_0_0/24 + rho**5/5760 + 5*rho**4/576 + rho**3*(17*z/216 + pi**2*z/12 + 7/108 + pi**2/12)/(6*(z + 1)) - rho**2*(-5*pi**2/72 + 5*zeta(3)/12)/2 + rho*(-11*z*zeta(3)/32 + pi**2*z/27 + 35*pi**4*z/864 - 11*zeta(3)/32 + 5*pi**2/144 + 35*pi**4/864)/(z + 1) - 307*pi**2*zeta(3)/864 - 19*zeta(5)/6 - 65*zeta(3)/648 + 53*pi**2/2916 + 365*pi**4/10368
```

### Af contribution

Consumed field: `/channels/q/I_coefficients`. External scalar multiplier: `1`; physical assembly sign: +1.

`eps^-6`:

```text
0
```

`eps^-5`:

```text
0
```

`eps^-4`:

```text
0
```

`eps^-3`:

```text
pi**2/432
```

`eps^-2`:

```text
-pi**2*H_z_m1/48 - H_z_m1_0_0/24 + rho**3/288 + pi**2*rho/96 - 5*zeta(3)/72 + pi**2/162
```

`eps^-1`:

```text
-3*pi**2*H_z_0_m1/32 - 3*H_z_0_m1_0_0/16 - H_z_m1*(-3*zeta(3)/8 + 13*pi**2/288) - 3*pi**2*H_z_m1_0/32 - 13*H_z_m1_0_0/144 - 3*H_z_m1_0_0_0/16 + 3*pi**2*H_z_m1_m1/16 + 3*H_z_m1_m1_0_0/8 + rho**4/192 - rho**3*(-z/18 - 5/144)/(6*z + 6) + pi**2*rho**2/48 + rho*(-3*z*zeta(3)/16 + 7*pi**2*z/288 - 3*zeta(3)/16 + pi**2/48)/(z + 1) - 83*zeta(3)/432 + 13*pi**2/972 + 11*pi**4/540
```

`eps^0`:

```text
-13*pi**2*H_z_0_0_m1/48 - 13*H_z_0_0_m1_0_0/24 - H_z_0_m1*(-41*zeta(3)/48 + 11*pi**2/64) - 25*pi**2*H_z_0_m1_0/144 - 11*H_z_0_m1_0_0/32 - 3*H_z_0_m1_0_0_0/16 + 41*pi**2*H_z_0_m1_m1/96 + 41*H_z_0_m1_m1_0_0/48 - H_z_m1*(-11*zeta(3)/16 + 31*pi**2/432 + 77*pi**4/1440) - H_z_m1_0*(-41*zeta(3)/48 + 11*pi**2/64) - H_z_m1_0_0*(31/216 + 43*pi**2/288) - 11*H_z_m1_0_0_0/32 - H_z_m1_0_0_0_0/6 + 41*pi**2*H_z_m1_0_m1/96 + 41*H_z_m1_0_m1_0_0/48 - H_z_m1_m1*(-11*pi**2/32 + 41*zeta(3)/24) + 41*pi**2*H_z_m1_m1_0/96 + 11*H_z_m1_m1_0_0/16 + 41*H_z_m1_m1_0_0_0/48 - 41*pi**2*H_z_m1_m1_m1/48 - 41*H_z_m1_m1_m1_0_0/24 + rho**5/5760 + 5*rho**4/576 + rho**3*(17*z/216 + pi**2*z/12 + 7/108 + pi**2/12)/(6*z + 6) - rho**2*(-5*pi**2/72 + 5*zeta(3)/12)/2 + rho*(-11*z*zeta(3)/32 + pi**2*z/27 + 35*pi**4*z/864 - 11*zeta(3)/32 + 5*pi**2/144 + 35*pi**4/864)/(z + 1) - 307*pi**2*zeta(3)/864 - 19*zeta(5)/6 - 65*zeta(3)/648 + 53*pi**2/2916 + 365*pi**4/10368
```

### Source-scoped symbol definitions

These definitions retain the original source notation; the normalized coefficient display above uses rho for its geometry logarithm.

```json
{
  "coefficient_notation": {
    "language": "Maple exact rational expressions",
    "argument": "z",
    "alphabet": [
      0,
      -1
    ],
    "Hlog(z,[a1,...,an])": "G(a1,...,an;z), kernels dt/(t-a1)",
    "Hlog(z,[-1])": "log(1+z)",
    "Hlog(z,[0])": "log(z)",
    "zeta[k]": "Riemann zeta(k)",
    "^": "power"
  }
}
```

Source dependency: [benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-04/derivation.md](sources/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-04/derivation.md).
Source dependency: [workflow/checks/shared-sign-production/publication-check.json](sources/workflow/checks/shared-sign-production/publication-check.json).
Source dependency: [benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-04/derivation.md](sources/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-04/derivation.md).
