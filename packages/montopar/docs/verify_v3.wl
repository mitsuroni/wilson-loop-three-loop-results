(* verify_v3.wl -- Use explicit FeynCalc Pair[Momentum[...],Momentum[...]] everywhere *)

$FeynCalcStartupMessages = False;
$FCCheckContext = False;
Needs["FeynCalc`"];

Print["===== VERIFICATION OF 4e-3g+3g GENERAL DERIVATION (v3) ====="];

(* Define shorthand for scalar products *)
sp[a_, b_] := Pair[Momentum[a], Momentum[b]];

(* TGV kinematic structures *)
GvA[nu1_, nu2_, nu5_] :=
  MT[nu1, nu2] * FV[l2 - l1, nu5] +
  MT[nu2, nu5] * FV[-l1 - 2*l2, nu1] +
  MT[nu5, nu1] * FV[2*l1 + l2, nu2];

GvB[nu3_, nu4_, nu6_] :=
  MT[nu3, nu4] * FV[l4 - l3, nu6] +
  MT[nu4, nu6] * FV[-l1 - l2 - l4, nu3] +
  MT[nu6, nu3] * FV[l1 + l2 + l3, nu4];

(* ============================================================ *)
Print[""];
Print["===== PART 1: General V_{2TGV} ====="];

(* FeynCalc direct *)
VfullFC = FV[vi, LorentzIndex[nu1]] * FV[vj, LorentzIndex[nu2]] *
          FV[vk, LorentzIndex[nu3]] * FV[vl, LorentzIndex[nu4]] *
          MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
          GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
          GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

Vfc = Contract[VfullFC] // ExpandScalarProduct;

(* Document's expression: G1 + G2 + G3 + G4a-d *)
(* Using sp[a,b] = Pair[Momentum[a], Momentum[b]] for consistency *)

G1 = sp[vi, vj] * sp[vk, vl] *
     (sp[l2, l4] - sp[l2, l3] - sp[l1, l4] + sp[l1, l3]);

G2 = sp[vi, vj] * (
  2*sp[vk, l1]*sp[vl, l2] - 2*sp[vk, l2]*sp[vl, l1] +
  (sp[vk, l2] - sp[vk, l1])*sp[vl, l3] +
  (sp[vl, l1] - sp[vl, l2])*sp[vk, l4]
);

G3 = 2*sp[vk, vl] * (
  (-sp[vi, l1] - 2*sp[vi, l2])*(sp[vj, l4] - sp[vj, l3]) +
  (2*sp[vj, l1] + sp[vj, l2])*(sp[vi, l4] - sp[vi, l3])
);

G4a = sp[vi, vk] * (2*sp[vj, l1] + sp[vj, l2]) *
      (sp[vl, l3] - sp[vl, l1] - sp[vl, l2]);
G4b = sp[vi, vl] * (2*sp[vj, l1] + sp[vj, l2]) *
      (sp[vk, l1] + sp[vk, l2] - sp[vk, l4]);
G4c = -sp[vj, vk] * (sp[vi, l1] + 2*sp[vi, l2]) *
      (sp[vl, l3] - sp[vl, l1] - sp[vl, l2]);
G4d = -sp[vj, vl] * (sp[vi, l1] + 2*sp[vi, l2]) *
      (sp[vk, l1] + sp[vk, l2] - sp[vk, l4]);

Vdoc = (G1 + G2 + G3 + G4a + G4b + G4c + G4d) // ExpandScalarProduct;

diffGen = Vfc - Vdoc // Simplify;
If[diffGen === 0,
  Print["✓ GENERAL EXPRESSION MATCHES"],
  Print["✗ GENERAL MISMATCH (non-zero after Simplify)"];
  Print["Difference has ", Length[diffGen // Expand], " terms after Expand"]
];

(* ============================================================ *)
Print[""];
Print["===== PART 2: S-channel (3,3,4,4) ====="];

ScalarProduct[v3, v3] = 0;
ScalarProduct[v4, v4] = 0;

(* FeynCalc direct S-channel *)
VfullS = FV[v3, LorentzIndex[nu1]] * FV[v3, LorentzIndex[nu2]] *
         FV[v4, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcS = Contract[VfullS] // ExpandScalarProduct;

(* Document general → S-channel substitution *)
VdocSubS = (Vdoc /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;

diffS = VfcS - VdocSubS // Simplify;
If[diffS === 0,
  Print["✓ S-CHANNEL: General→S substitution matches FeynCalc"],
  Print["✗ S-CHANNEL: General→S substitution MISMATCH"];
  Print["Difference: ", diffS]
];

(* Expected S-channel: z (v3.(l1-l2)) (v4.(l3-l4)) *)
expectedS = sp[v3, v4] * (sp[v3, l1] - sp[v3, l2]) * (sp[v4, l3] - sp[v4, l4]);

diffSexp = VfcS - expectedS // Simplify;
If[diffSexp === 0,
  Print["✓ S-CHANNEL: FeynCalc = z(v3.(l1-l2))(v4.(l3-l4))"],
  Print["✗ S-CHANNEL: FeynCalc ≠ expected"];
  Print["FeynCalc: ", VfcS];
  Print["Expected: ", expectedS];
  Print["Diff: ", diffSexp]
];

(* ============================================================ *)
Print[""];
Print["===== PART 3: S-channel degeneration ====="];

G1s = (G1 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // Simplify;
G2s = (G2 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // Simplify;
G3s = (G3 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // Simplify;
G4sum = (G4a + G4b + G4c + G4d /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct // Simplify;

Print["G1 → ", If[G1s===0, "✓ 0", "✗ "<>ToString[G1s]]];
Print["G2 → ", If[G2s===0, "✓ 0", "✗ "<>ToString[G2s]]];
Print["G3 → ", If[G3s===0, "✓ 0", "✗ "<>ToString[G3s]]];
diffG4 = G4sum - expectedS // Simplify;
Print["G4a-d sum vs expected: ", If[diffG4===0, "✓ MATCHES", "✗ DIFF: "<>ToString[diffG4]]];

(* ============================================================ *)
Print[""];
Print["===== PART 4: T-channel (3,4,3,4) ====="];

VfullT = FV[v3, LorentzIndex[nu1]] * FV[v4, LorentzIndex[nu2]] *
         FV[v3, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcT = Contract[VfullT] // ExpandScalarProduct;

VdocSubT = (Vdoc /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // ExpandScalarProduct;

diffT = VfcT - VdocSubT // Simplify;
If[diffT === 0,
  Print["✓ T-CHANNEL: General→T substitution matches FeynCalc"],
  Print["✗ T-CHANNEL: General→T substitution MISMATCH"];
  Print["FeynCalc T: ", VfcT];
  Print[""];
  Print["Document T: ", VdocSubT];
  Print[""];
  Print["Difference: ", diffT]
];

(* ============================================================ *)
Print[""];
Print["===== PART 5: Verify G1+G2+G3+G4 against FeynCalc for T-channel ====="];

(* Check each group separately for T-channel *)
G1t = (G1 /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // Simplify;
G2t = (G2 /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // Simplify;
G3t = (G3 /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // Simplify;
G4t = (G4a + G4b + G4c + G4d /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // Simplify;

Print["T-channel G1: ", G1t];
Print["T-channel G2: ", G2t];
Print["T-channel G3: ", G3t];
Print["T-channel G4: ", G4t];

Print[""];
Print["===== VERIFICATION COMPLETE ====="];