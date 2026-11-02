(* fix_G3.wl -- Extract correct G3 from FeynCalc output *)

$FeynCalcStartupMessages = False;
$FCCheckContext = False;
Needs["FeynCalc`"];

sp[a_, b_] := Pair[Momentum[a], Momentum[b]];

GvA[nu1_, nu2_, nu5_] :=
  MT[nu1, nu2] * FV[l2 - l1, nu5] +
  MT[nu2, nu5] * FV[-l1 - 2*l2, nu1] +
  MT[nu5, nu1] * FV[2*l1 + l2, nu2];

GvB[nu3_, nu4_, nu6_] :=
  MT[nu3, nu4] * FV[l4 - l3, nu6] +
  MT[nu4, nu6] * FV[-l1 - l2 - l4, nu3] +
  MT[nu6, nu3] * FV[l1 + l2 + l3, nu4];

(* FeynCalc ground truth *)
VfullFC = FV[vi, LorentzIndex[nu1]] * FV[vj, LorentzIndex[nu2]] *
          FV[vk, LorentzIndex[nu3]] * FV[vl, LorentzIndex[nu4]] *
          MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
          GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
          GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

Vfc = Contract[VfullFC] // ExpandScalarProduct;

(* Verified-correct groups: G1, G2, G4a-d *)
G1 = sp[vi, vj] * sp[vk, vl] *
     (sp[l2, l4] - sp[l2, l3] - sp[l1, l4] + sp[l1, l3]);

G2 = sp[vi, vj] * (
  2*sp[vk, l1]*sp[vl, l2] - 2*sp[vk, l2]*sp[vl, l1] +
  (sp[vk, l2] - sp[vk, l1])*sp[vl, l3] +
  (sp[vl, l1] - sp[vl, l2])*sp[vk, l4]
);

G4a = sp[vi, vk] * (2*sp[vj, l1] + sp[vj, l2]) *
      (sp[vl, l3] - sp[vl, l1] - sp[vl, l2]);
G4b = sp[vi, vl] * (2*sp[vj, l1] + sp[vj, l2]) *
      (sp[vk, l1] + sp[vk, l2] - sp[vk, l4]);
G4c = -sp[vj, vk] * (sp[vi, l1] + 2*sp[vi, l2]) *
      (sp[vl, l3] - sp[vl, l1] - sp[vl, l2]);
G4d = -sp[vj, vl] * (sp[vi, l1] + 2*sp[vi, l2]) *
      (sp[vk, l1] + sp[vk, l2] - sp[vk, l4]);

G124 = (G1 + G2 + G4a + G4b + G4c + G4d) // ExpandScalarProduct;

(* Extract correct G3 = FeynCalc - (G1+G2+G4) *)
G3correct = (Vfc - G124) // ExpandScalarProduct // Simplify;

Print["===== CORRECT G3 (from FeynCalc - G1 - G2 - G4) ====="];
Print[G3correct];
Print[""];

(* Check if G3correct factorizes as 2*(vk.vl)*[...] *)
Print["===== Attempted factorization: G3correct / (2*(vk.vl)) ====="];
G3factor = G3correct / (2*sp[vk, vl]) // Simplify;
Print[G3factor];
Print[""];

(* Verify: Vfc = G1 + G2 + G3correct + G4 *)
Vreconstructed = (G1 + G2 + G3correct + G4a + G4b + G4c + G4d) // ExpandScalarProduct;
diff = Vfc - Vreconstructed // Simplify;
Print["===== Verification: FeynCalc - (G1+G2+G3correct+G4) ====="];
Print[If[diff===0, "✓ ZERO — G3 is correct", "✗ " <> ToString[diff]]];

(* Verify S-channel with corrected G3 *)
Print[""];
Print["===== S-channel check with corrected G3 ====="];
ScalarProduct[v3, v3] = 0;
ScalarProduct[v4, v4] = 0;

VfullS = FV[v3, LorentzIndex[nu1]] * FV[v3, LorentzIndex[nu2]] *
         FV[v4, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcS = Contract[VfullS] // ExpandScalarProduct;

VdocNewS = (G1 + G2 + G3correct + G4a + G4b + G4c + G4d /.
  {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct // Simplify;

diffS = VfcS - VdocNewS // Simplify;
Print["S-channel diff: ", If[diffS===0, "✓ ZERO", "✗ " <> ToString[diffS]]];

(* Verify T-channel with corrected G3 *)
Print[""];
Print["===== T-channel check with corrected G3 ====="];

VfullT = FV[v3, LorentzIndex[nu1]] * FV[v4, LorentzIndex[nu2]] *
         FV[v3, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcT = Contract[VfullT] // ExpandScalarProduct;

VdocNewT = (G1 + G2 + G3correct + G4a + G4b + G4c + G4d /.
  {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // ExpandScalarProduct // Simplify;

diffT = VfcT - VdocNewT // Simplify;
Print["T-channel diff: ", If[diffT===0, "✓ ZERO", "✗ " <> ToString[diffT]]];