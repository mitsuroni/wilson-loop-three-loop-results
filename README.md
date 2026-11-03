# Three-loop Wilson-loop results in QCD

Analytic results for the connected lightlike quadrilateral Wilson-loop exponent in massless QCD, through three loops and through the finite term.

**[Four-color analytic formulas](results/current/four-color-results.md)** · **[Current exact data](benchmarks/current/README.md)** · **[Integral derivations](integrals/README.md)** · **[Reproduce the calculation](docs/reproduction.md)**

With $T=T_Fn_f$ and $\lambda=\alpha_s^{\overline{\mathrm{MS}}}/\pi$, the three-loop result is organized as

$$F_3=C_F\left[C_A^2F_{3,AA}+C_ATF_{3,Af}+C_FTF_{3,Ff}+T^2F_{3,ff}\right].$$

The [result article](results/current/four-color-results.md) displays every Laurent coefficient of the scale-stripped sum and bare exponent, the coupling-renormalized poles, and all four finite functions $F_{3,c}(L,z)$ and $f_{3,c}(z)=F_{3,c}(0,z)$. It defines the functions, normalization and scale convention alongside the formulas.

## Start with your research question

| You want to… | Open |
|---|---|
| Read the final analytic formulas | [Four-color Laurent coefficients and finite terms](results/current/four-color-results.md) |
| Load the complete exact result | [Three-loop JSON](benchmarks/current/results/3loop.json); [one-loop](benchmarks/current/results/1loop.json) and [two-loop](benchmarks/current/results/2loop.json) inputs |
| Find an original graph or C/J integral | [Current163-source catalogue](benchmarks/current/catalog.json), [prepared inputs](benchmarks/current/prepared/), [original input catalogue](data/README.md) |
| Follow an individual integral calculation | [Integral reading index](integrals/README.md) |
| Run source assembly, counterterms and RG checks | [Current calculation guide](docs/current/README.md), [reproduction commands](docs/reproduction.md) |
| Inspect checks and mathematical provenance | [Validation index](validation/README.md), [single convention](workflow/current-conventions.md) |

## Current data and derivations

The current catalogue covers **163 source orbits**:2 at one loop,12 at two loops and149 at three loops. Each selected integral is connected to its graph and native C/J definition. Four-leg sources are assembled through the joint-object ledger.

The continuous reading layer covers **2 one-loop,12 two-loop and39 three-loop6e integrals**. Each entry connects its actual initial integral to evaluated subintegrals and the requested Laurent coefficients. The [15 insertion](evaluations/three-loop-insertions/README.md) and [44 tree-three-gluon](evaluations/three-loop-tree5e/README.md) calculation packets are also available. Their original methods, code, raw outputs and validation records accompany the results.

The [validation index](validation/README.md) connects the complete four-color assembly and RG checks, the53-entry public replay, and the individual mathematical reading audits.

## Run the current calculation

From the repository root:

```text
python -m pip install -r docs/current/requirements.txt
python -m workflow.current_calculation check
python -m workflow.current_calculation result 3
python -m workflow.current_calculation replay --output-dir replay-output/current-new
```

Choose a new output directory for each replay. The commands use the source and exact-data snapshot included in this repository. See [reproduction](docs/reproduction.md) for individual integral batches and formula generation.

## Repository map

| Directory | Contents |
|---|---|
| [results/current](results/current/README.md) | Main result article and exact-data entrypoints |
| [benchmarks/current](benchmarks/current/README.md) | Current source catalogue, prepared C/J and one-to-three-loop result data |
| [integrals](integrals/README.md) | Individual derivations, shared mathematics, code and data indices |
| [evaluations](docs/README.md#calculation-packages) | Reproducible family calculation packets |
| [packages](packages/README.md), [data](data/README.md) | Diagram generator, Montopar and original graph/C/J snapshots |
| [validation](validation/README.md) | Calculation checks, exact formula checks and public-copy execution |
| [docs](docs/README.md), [workflow](workflow/current-conventions.md) | Reading guide, conventions and calculation interface |

## Historical archive

[Dated releases and research records](docs/archive/README.md) preserve the evolution of the calculation and its original evidence. Use the current result, data and calculation links above for present work.

Research materials prepared by Dengyao. Source provenance and original third-party notices accompany their respective files.
