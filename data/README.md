# Graph definitions and original integral inputs

Start with the fixed [2026-11-02 input catalogue](2026-11-02/README.md).
It contains all **163 source orbits** through three loops: 2 at one loop,
12 at two loops, and 149 at three loops. A source orbit is a canonical
representative in the adopted catalogue; an effective self-energy or vertex
operator is not counted as one newly generated microscopic QCD diagram.

- [Machine catalogue](2026-11-02/catalog.json): stable IDs, families, graph
  definitions, C/J references, normalization notes and four-leg grouping.
- [Coverage and explicit gaps](2026-11-02/coverage.json).
- [Byte-preserved source manifest](2026-11-02/source-manifest.json).
- [All publication file hashes](manifest.json).
- [Generator package](../packages/qgraf-pipeline/README.md) and
  [Montopar package](../packages/montopar/README.md).

Every integral has its own `input.json` index and English `README.md` below
`2026-11-02/integrals/<N>loop/<orbit_id>/`. The indexes point into complete,
unmodified source files; they do not replace a large polynomial, ordered
domain or block sum with a display summary. JSON pointer field hashes use
canonical UTF-8 JSON, with sorted keys, compact separators and no ASCII
escaping. File hashes always apply to the actual original bytes.

This is the first **upstream-data** publication: graph -> original C/J ->
available prepared inputs. It exports no new Laurent results and executes
no integral calculation. Previously calculated coefficients and the later
four-color acceptance have their own versioned result records; absence from
this input-only batch does not mean a coefficient is mathematically unknown.

For an offline integrity and pointer/navigation audit, run from the repository
root with the Python standard library:

```text
python tools/check_upstream.py --output-dir upstream-audit
```

Use a fresh audit output directory. This checks included files and field
bindings, not graph completeness or the correctness of an evaluated integral.
Actual source generation, symbolic integration and function-value evaluation
are separate reproduction levels.
