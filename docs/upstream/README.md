# From diagrams to parameter-integral inputs

This is the first curated part of the research repository. It follows the scientific data flow without requiring the full Wilson development environment.

1. Read the [diagram-generation guide](generation.md) and the [generator source](../../packages/qgraf-pipeline/README.md).
2. Open the [163-integral English index](catalogue.md) or the [frozen input catalogue](../../data/README.md), select an orbit, and follow its original definition.
3. Read the [Montopar guide](parameterization.md), then inspect the original C/J record and its variables, measure, bounds and representative.
4. Check the input's source hashes and normalization notes before comparing it with any evaluated integral.

The catalogue includes microscopic or effective-subgraph orbit records as identified by their source. An orbit count is not a count of distinct microscopic QCD diagrams. The four-leg raw graph sources and their later joint integrals remain separately identified.

## What this stage produces

Montopar outputs a parameter-integral definition. In the raw convention a contribution has the form:

```math
\mathcal A_\rho=C_{\rho,\mathrm{raw}}J_{\rho,\mathrm{raw}}.
```

An optional IBP decomposition instead carries a sum of terms, each with its own coefficient, integrand and domain. It does not supply an evaluated Laurent series. The next curation stage will add a continuous English evaluation article, the actual scripts/data, and checks for each independent integral or joint object.

See [publication status](../publication-status.md) and [reproduction levels](../reproduction.md).
