# Preserved low-loop calculation context

This shared directory keeps the byte-exact source/code/data context used by the previously published low-loop evaluators. Current reader articles are under integrals/one-loop and integrals/two-loop in the main project; older Markdown inside this runtime is archival provenance, not the new reader layer.

The relative subtrees retain their original context so that programs and input hashes do not need scientific patches. The original one-loop program calculates its own source geometry, two parameter integrals, Gamma/M expansion and checks. The two-loop family engine expands the evaluated Gamma/M mother formulas and adopts previously evaluated TGV04 HLOG/compact jets; its separate controls check original ordered kernels, Pi1 tensor/corners and real moment integrals. These roles are not12 new cold CAS integrations.

`source-manifest.json` records exact byte provenance for this shared runtime. Per-entry code/reproduce.py invokes the actual same-loop batch in a new output directory. Archived Maple/Wolfram scripts are evidence of original operations and may contain old output paths; they are not automatically launched. Native-language files are preserved, not created merely to fill a three-language template.

The main project's workflow/current-conventions.md remains the unique physical interface. This directory introduces no alternative source, colour or normalization convention and does not change accepted results. Runtime checks and the independent review of the new mathematical prose are reported separately.
