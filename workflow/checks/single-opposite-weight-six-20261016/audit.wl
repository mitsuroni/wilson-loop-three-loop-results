Module[{d, key, vec, all, coords, dense, cols, es, checks, witness, pair, r0, r1},
 d = Import["D:/wilson-agent/workflow/checks/single-opposite-weight-six-20261016/reduction.json", "RawJSON"];
 key[t_] := ToString[{t["word"], t["zeta_powers"]}, InputForm];
 vec[ts_] := Association[(key[#] -> ToExpression[#["coefficient"]]) & /@ ts];
 all = Join[(#["hlog_terms"] & /@ d["options"]), Lookup[d["moments"], {"200", "110", "020", "011"}]];
 coords = Union[Flatten[Keys[vec[#]] & /@ all]];
 dense[ts_] := Lookup[vec[ts], coords, 0];
 cols = Transpose[dense[#["hlog_terms"]] & /@ d["options"]];
 es = Transpose[dense /@ Lookup[d["moments"], {"200", "110", "020", "011"}]];
 r0 = MatrixRank[cols]; r1 = MatrixRank[Join[cols, es, 2]];
 checks = {
   r0 == 12, r1 == 13,
   MatrixRank[Join[cols, es[[All, {1}]], 2]] == 13,
   MatrixRank[Join[cols, es[[All, {2}]], 2]] == 13,
   MatrixRank[Join[cols, es[[All, {3}]], 2]] == 12,
   MatrixRank[Join[cols, es[[All, {4}]], 2]] == 13,
   MatrixRank[Join[cols, es[[All, {2}]] - es[[All, {1}]]/2, 2]] == 12,
   MatrixRank[Join[cols, es[[All, {4}]] - es[[All, {1}]]/2, 2]] == 12,
   MatrixRank[Join[cols, es[[All, {3}]] - es[[All, {4}]], 2]] == 13,
   MatrixRank[Join[cols, Transpose[{dense[d["I"][[3]]]}], 2]] == 13
 };
 witness = d["period_free_obstruction_E200"]["witness"];
 pair[ts_] := Total[(ToExpression[#["coefficient"]] Lookup[vec[ts], key[#["coordinate"]], 0]) & /@ witness];
 checks = Join[checks, {And @@ (pair[#["hlog_terms"]] == 0 & /@ d["options"]),
   (pair /@ Lookup[d["moments"], {"200", "110", "020", "011"}]) == {-2,-1,0,-1}}];
 Export["D:/wilson-agent/workflow/checks/single-opposite-weight-six-20261016/wolfbook-check.json",
   <|"passed" -> And @@ checks, "checks" -> Length[checks], "base_rank" -> r0,
     "joint_rank" -> r1, "period_free_values" -> (pair /@ Lookup[d["moments"], {"200", "110", "020", "011"}]),
     "formal" -> "not_run"|>, "RawJSON"];
 And @@ checks
]
