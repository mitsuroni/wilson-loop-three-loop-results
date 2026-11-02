(* Color trace template for montopar *)
(* INDICES = comma-separated SUNT color indices *)
(* TGV_FACTOR = structure constant factor for TGV vertices, or 1 *)

<< FeynCalc`

result = SUNSimplify[SUNT[{{INDICES}}] * {{TGV_FACTOR}}, SUNNToCACF -> True];
Print[InputForm[result]];
