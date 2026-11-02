(* verify_merged_form.wl -- Verify merged G3+G4 organization *)

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

VfullFC = FV[vi, LorentzIndex[nu1]] * FV[vj, LorentzIndex[nu2]] *
          FV[vk, LorentzIndex[nu3]] * FV[vl, LorentzIndex[nu4]] *
          MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
          GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
          GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

Vfc = Contract[VfullFC] // ExpandScalarProduct;

(* MERGED FORM: G3 cross terms absorbed into G4 *)
(* G1, G2 unchanged *)
G1 = sp[vi, vj] * sp[vk, vl] *
     (sp[l2, l4] - sp[l2, l3] - sp[l1, l4] + sp[l1, l3]);

G2 = sp[vi, vj] * (
  2*sp[vk, l1]*sp[vl, l2] - 2*sp[vk, l2]*sp[vl, l1] +
  (sp[vk, l2] - sp[vk, l1])*sp[vl, l3] +
  (sp[vl, l1] - sp[vl, l2])*sp[vk, l4]
);

(* G3 reduced: only (vi.vj) and (vk.vl) terms *)
G3_reduced = sp[vi, vj] * 4*(sp[vl, l1]*sp[vk, l2] - sp[vk, l1]*sp[vl, l2]) +
             sp[vk, vl] * ((sp[vi, l1] + 2*sp[vi, l2])*(sp[vj, l3] - sp[vj, l4]) -
                           (2*sp[vj, l1] + sp[vj, l2])*(sp[vi, l3] - sp[vi, l4]));

(* G4 merged: original G4a-d + G3 cross terms *)
G4a_new = sp[vi, vk] * (2*sp[vj, l1] + sp[vj, l2]) *
          (sp[vl, l1] + sp[vl, l2] + sp[vl, l3]);
G4b_new = sp[vi, vl] * (2*sp[vj, l1] + sp[vj, l2]) *
          (-sp[vk, l1] - sp[vk, l2] - sp[vk, l4]);
G4c_new = -sp[vj, vk] * (sp[vi, l1] + 2*sp[vi, l2]) *
          (sp[vl, l1] + sp[vl, l2] + sp[vl, l3]);
G4d_new = -sp[vj, vl] * (sp[vi, l1] + 2*sp[vi, l2]) *
          (-sp[vk, l1] - sp[vk, l2] - sp[vk, l4]);

Vmerged = (G1 + G2 + G3_reduced + G4a_new + G4b_new + G4c_new + G4d_new) // ExpandScalarProduct;

diff = Vfc - Vmerged // Simplify;
Print["===== Merged form verification ====="];
Print["FeynCalc = G1 + G2 + G3_reduced + G4_new ?"];
Print[If[diff===0, "✓ YES", "✗ NO, diff: " <> ToString[diff]]];

(* S and T channel checks *)
Print[""];
Print["===== S-channel ====="];
ScalarProduct[v3, v3] = 0;
ScalarProduct[v4, v4] = 0;

VfullS = FV[v3, LorentzIndex[nu1]] * FV[v3, LorentzIndex[nu2]] *
         FV[v4, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];
VfcS = Contract[VfullS] // ExpandScalarProduct;
VmergedS = (Vmerged /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct // Simplify;
Print[If[VfcS-VmergedS//Simplify===0, "✓ MATCHES", "✗ MISMATCH"]];

Print[""];
Print["===== T-channel ====="];
VfullT = FV[v3, LorentzIndex[nu1]] * FV[v4, LorentzIndex[nu2]] *
         FV[v3, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];
VfcT = Contract[VfullT] // ExpandScalarProduct;
VmergedT = (Vmerged /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // ExpandScalarProduct // Simplify;
Print[If[VfcT-VmergedT//Simplify===0, "✓ MATCHES", "✗ MISMATCH"]];

Print[""];
Print["===== CORRECT V_{2TGV} (MERGED FORM) ====="];
Print[""];
Print["G1 = (vi.vj)(vk.vl) * (l2-l1).(l4-l3)"];
Print[""];
Print["G2 = (vi.vj) * ["];
Print["  2(vk.l1)(vl.l2) - 2(vk.l2)(vl.l1)"];
Print["  + (vk.l2 - vk.l1)(vl.l3)"];
Print["  + (vl.l1 - vl.l2)(vk.l4)"];
Print["]"];
Print[""];
Print["G3 = (vi.vj) * 4[(vl.l1)(vk.l2) - (vk.l1)(vl.l2)]"];
Print["  + (vk.vl) * [(vi.l1+2vi.l2)(vj.l3-vj.l4)"];
Print["                - (2vj.l1+vj.l2)(vi.l3-vi.l4)]"];
Print[""];
Print["G4a = (vi.vk) * (2vj.l1+vj.l2) * (vl.l1 + vl.l2 + vl.l3)"];
Print["G4b = (vi.vl) * (2vj.l1+vj.l2) * (-vk.l1 - vk.l2 - vk.l4)"];
Print["G4c = -(vj.vk) * (vi.l1+2vi.l2) * (vl.l1 + vl.l2 + vl.l3)"];
Print["G4d = -(vj.vl) * (vi.l1+2vi.l2) * (-vk.l1 - vk.l2 - vk.l4)"];