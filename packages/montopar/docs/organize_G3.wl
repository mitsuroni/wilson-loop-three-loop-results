(* organize_G3.wl -- Organize correct G3 into readable form *)

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

(* Verified groups *)
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

G3correct = (Vfc - (G1 + G2 + G4a + G4b + G4c + G4d)) // ExpandScalarProduct // Simplify;

(* Try to factor G3correct by grouping *)
(* Collect terms by velocity-velocity prefactor *)

Print["===== CORRECT V_{2TGV} (FeynCalc ground truth) ====="];
Print[""];
Print["G1 = (vi.vj)(vk.vl) (l2-l1).(l4-l3)"];
Print[""];
Print["G2 = (vi.vj) [ 2(vk.l1)(vl.l2) - 2(vk.l2)(vl.l1) + (vk.l2-vk.l1)(vl.l3) + (vl.l1-vl.l2)(vk.l4) ]"];
Print[""];
Print["G3 = CORRECTED (was wrong in old doc):"];
Print[""];

(* G3 has terms with all 6 v.v prefactors. Let's see if it simplifies *)
(* Try factoring each v.v type *)

(* Collect terms with vi.vj *)
G3_vivj = Collect[G3correct, sp[vi, vj]] // Simplify;
Print["G3 terms with (vi.vj) prefactor:"];
Print[G3_vivj];
Print[""];

(* Collect terms with vk.vl *)
G3_vkvl = Collect[G3correct, sp[vk, vl]] // Simplify;
Print["G3 terms with (vk.vl) prefactor:"];
Print[G3_vkvl];
Print[""];

(* Check if G3 = 2*(vk.vl)*[...] + (vi.vj)*[...] + cross terms *)
(* Try factorization *)
G3a = Coefficient[G3correct, sp[vi, vj]];
G3b = Coefficient[G3correct, sp[vk, vl]];
G3c = Coefficient[G3correct, sp[vi, vk]];
G3d = Coefficient[G3correct, sp[vi, vl]];
G3e = Coefficient[G3correct, sp[vj, vk]];
G3f = Coefficient[G3correct, sp[vj, vl]];

Print["===== G3 factorization by v.v prefactor ====="];
Print["(vi.vj) coefficient: ", G3a];
Print["(vk.vl) coefficient: ", G3b];
Print["(vi.vk) coefficient: ", G3c];
Print["(vi.vl) coefficient: ", G3d];
Print["(vj.vk) coefficient: ", G3e];
Print["(vj.vl) coefficient: ", G3f];