(* verify_O_v2.wl -- Fixed: replace at Momentum level *)

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

(* ============================================================ *)
(* Reformulate L-to-O conversion at Momentum level *)
(* ============================================================ *)

(* FeynCalc V_{2TGV} *)
VfullFC = FV[vi, LorentzIndex[nu1]] * FV[vj, LorentzIndex[nu2]] *
          FV[vk, LorentzIndex[nu3]] * FV[vl, LorentzIndex[nu4]] *
          MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
          GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
          GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

Vfc = Contract[VfullFC] // ExpandScalarProduct;

(* Replace Momentum[l_i] -> I * Momentum[dX_i] *)
(* This automatically handles Pair[Momentum[l_i], ...] in all positions *)
OfromFC = Vfc /. {
  Momentum[l1] -> I * Momentum[dX1],
  Momentum[l2] -> I * Momentum[dX2],
  Momentum[l3] -> I * Momentum[dX3],
  Momentum[l4] -> I * Momentum[dX4]
};
OfromFC = OfromFC // ExpandScalarProduct // Simplify;

Print["===== PART 1: General O operator ====="];

(* Document's O (same as before, verified) *)
OI = -sp[vi, vj] * sp[vk, vl] *
     (sp[dX2, dX4] - sp[dX2, dX3] - sp[dX1, dX4] + sp[dX1, dX3]);

OIIabG2 = -sp[vi, vj] * (
  2*sp[vk, dX1]*sp[vl, dX2] - 2*sp[vk, dX2]*sp[vl, dX1] +
  (sp[vk, dX2] - sp[vk, dX1])*sp[vl, dX3] +
  (sp[vl, dX1] - sp[vl, dX2])*sp[vk, dX4]
);

OIIabG3 = -4*sp[vi, vj] * (sp[vl, dX1]*sp[vk, dX2] - sp[vk, dX1]*sp[vl, dX2]);

OIIcd = -sp[vk, vl] * (
  (sp[vi, dX1] + 2*sp[vi, dX2])*(sp[vj, dX3] - sp[vj, dX4]) -
  (2*sp[vj, dX1] + sp[vj, dX2])*(sp[vi, dX3] - sp[vi, dX4])
);

OIIIg3 = -2*sp[vi, vk]*(2*sp[vj, dX1]+sp[vj, dX2])*(sp[vl, dX1]+sp[vl, dX2]) +
          2*sp[vi, vl]*(2*sp[vj, dX1]+sp[vj, dX2])*(sp[vk, dX1]+sp[vk, dX2]) +
          2*sp[vj, vk]*(sp[vi, dX1]+2*sp[vi, dX2])*(sp[vl, dX1]+sp[vl, dX2]) -
          2*sp[vj, vl]*(sp[vi, dX1]+2*sp[vi, dX2])*(sp[vk, dX1]+sp[vk, dX2]);

OIIIg4 = -sp[vi, vk]*(2*sp[vj, dX1]+sp[vj, dX2])*(sp[vl, dX3]-sp[vl, dX1]-sp[vl, dX2]) -
          sp[vi, vl]*(2*sp[vj, dX1]+sp[vj, dX2])*(sp[vk, dX1]+sp[vk, dX2]-sp[vk, dX4]) +
          sp[vj, vk]*(sp[vi, dX1]+2*sp[vi, dX2])*(sp[vl, dX3]-sp[vl, dX1]-sp[vl, dX2]) +
          sp[vj, vl]*(sp[vi, dX1]+2*sp[vi, dX2])*(sp[vk, dX1]+sp[vk, dX2]-sp[vk, dX4]);

Odoc = (OI + OIIabG2 + OIIabG3 + OIIcd + OIIIg3 + OIIIg4) // ExpandScalarProduct // Simplify;

diff = OfromFC - Odoc // Simplify;
If[Simplify[diff] == 0,
  Print["✓ GENERAL O OPERATOR MATCHES"],
  Print["✗ MISMATCH. Difference:"];
  Print[diff]
];

(* ============================================================ *)
Print[""];
Print["===== PART 2: S-channel O reduction ====="];

ScalarProduct[v3, v3] = 0;
ScalarProduct[v4, v4] = 0;

VfullS = FV[v3, LorentzIndex[nu1]] * FV[v3, LorentzIndex[nu2]] *
         FV[v4, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcS = Contract[VfullS] // ExpandScalarProduct;

OfromFCS = (VfcS /. {
  Momentum[l1] -> I * Momentum[dX1],
  Momentum[l2] -> I * Momentum[dX2],
  Momentum[l3] -> I * Momentum[dX3],
  Momentum[l4] -> I * Momentum[dX4]
}) // ExpandScalarProduct // Simplify;

OdocS = (Odoc /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct // Simplify;

expectedOS = -sp[v3, v4] * (sp[v3, dX1] - sp[v3, dX2]) * (sp[v4, dX3] - sp[v4, dX4]);

Print["O_S (FeynCalc):  ", OfromFCS];
Print["O_S (Document):  ", OdocS];
Print["O_S (Expected):  ", expectedOS];
Print[""];
Print["FeynCalc vs Expected: ", If[Simplify[OfromFCS - expectedOS] == 0, "✓", "✗"]];
Print["Document vs Expected: ", If[Simplify[OdocS - expectedOS] == 0, "✓", "✗"]];
Print["FeynCalc vs Document: ", If[Simplify[OfromFCS - OdocS] == 0, "✓", "✗"]];

(* ============================================================ *)
Print[""];
Print["===== PART 3: T-channel O reduction ====="];

VfullT = FV[v3, LorentzIndex[nu1]] * FV[v4, LorentzIndex[nu2]] *
         FV[v3, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcT = Contract[VfullT] // ExpandScalarProduct;

OfromFCT = (VfcT /. {
  Momentum[l1] -> I * Momentum[dX1],
  Momentum[l2] -> I * Momentum[dX2],
  Momentum[l3] -> I * Momentum[dX3],
  Momentum[l4] -> I * Momentum[dX4]
}) // ExpandScalarProduct // Simplify;

OdocT = (Odoc /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // ExpandScalarProduct // Simplify;

diffT = OfromFCT - OdocT // Simplify;
Print["O_T FeynCalc vs Document: ", If[Simplify[diffT] == 0, "✓ MATCHES", "✗ DIFF: " <> ToString[diffT]]];

Print[""];
Print["===== VERIFICATION COMPLETE ====="];