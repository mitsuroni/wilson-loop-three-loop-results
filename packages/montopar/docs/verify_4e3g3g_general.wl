(* verify_4e3g3g_general.wl -- Comprehensive verification of V_{2TGV} *)
(* Compares document's 9-term expression with FeynCalc raw output *)
(* Also verifies S-channel and T-channel specializations *)

$FeynCalcStartupMessages = False;
$FCCheckContext = False;
Needs["FeynCalc`"];

Print["===== VERIFICATION OF 4e-3g+3g GENERAL DERIVATION ====="];
Print[""];

(* ============================================================ *)
Print["===== PART 1: General V_{2TGV} Lorentz contraction ====="];
Print[""];

(* TGV A: G^A_{nu1,nu2,nu5} with k1=-l1, k2=-l2, k3=l1+l2 *)
GvA[nu1_, nu2_, nu5_] :=
  MT[nu1, nu2] * FV[l2 - l1, nu5] +
  MT[nu2, nu5] * FV[-l1 - 2*l2, nu1] +
  MT[nu5, nu1] * FV[2*l1 + l2, nu2];

(* TGV B: G^B_{nu3,nu4,nu6} with k1=-l3, k2=-l4, k3=l1+l2 *)
GvB[nu3_, nu4_, nu6_] :=
  MT[nu3, nu4] * FV[l4 - l3, nu6] +
  MT[nu4, nu6] * FV[-l1 - l2 - l4, nu3] +
  MT[nu6, nu3] * FV[l1 + l2 + l3, nu4];

(* FeynCalc: general contraction with symbolic vi, vj, vk, vl *)
VfullFC = FV[vi, LorentzIndex[nu1]] * FV[vj, LorentzIndex[nu2]] *
          FV[vk, LorentzIndex[nu3]] * FV[vl, LorentzIndex[nu4]] *
          MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
          GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
          GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VcontractedFC = Contract[VfullFC];
VexpandedFC = ExpandScalarProduct[VcontractedFC];

Print["FeynCalc raw output (general vi,vj,vk,vl):"];
Print[VexpandedFC];
Print[""];

(* Document's 9-term expression, organized into 4 groups *)
(* Using symbolic vi, vj, vk, vl for the velocity vectors *)
(* G1: (vi.vj)(vk.vl) (l2-l1).(l4-l3) *)
G1 = Pair[vi, vj] * Pair[vk, vl] *
     (Pair[l2, l4] - Pair[l2, l3] - Pair[l1, l4] + Pair[l1, l3]);

(* G2: (vi.vj)[...] *)
G2 = Pair[vi, vj] * (
  2*Pair[vk, l1]*Pair[vl, l2] - 2*Pair[vk, l2]*Pair[vl, l1] +
  (Pair[vk, l2] - Pair[vk, l1])*Pair[vl, l3] +
  (Pair[vl, l1] - Pair[vl, l2])*Pair[vk, l4]
);

(* G3: 2(vk.vl)[...] *)
G3 = 2*Pair[vk, vl] * (
  (-Pair[vi, l1] - 2*Pair[vi, l2])*(Pair[vj, l4] - Pair[vj, l3]) +
  (2*Pair[vj, l1] + Pair[vj, l2])*(Pair[vi, l4] - Pair[vi, l3])
);

(* G4a-d *)
G4a = Pair[vi, vk] * (2*Pair[vj, l1] + Pair[vj, l2]) *
      (Pair[vl, l3] - Pair[vl, l1] - Pair[vl, l2]);
G4b = Pair[vi, vl] * (2*Pair[vj, l1] + Pair[vj, l2]) *
      (Pair[vk, l1] + Pair[vk, l2] - Pair[vk, l4]);
G4c = -Pair[vj, vk] * (Pair[vi, l1] + 2*Pair[vi, l2]) *
      (Pair[vl, l3] - Pair[vl, l1] - Pair[vl, l2]);
G4d = -Pair[vj, vl] * (Pair[vi, l1] + 2*Pair[vi, l2]) *
      (Pair[vk, l1] + Pair[vk, l2] - Pair[vk, l4]);

Vdoc = G1 + G2 + G3 + G4a + G4b + G4c + G4d;
VdocExpanded = ExpandScalarProduct[Vdoc];

Print["Document's expression (general vi,vj,vk,vl):"];
Print[VdocExpanded];
Print[""];

(* Compare *)
diff = ExpandScalarProduct[VexpandedFC - VdocExpanded];
Print["Difference (FeynCalc - Document):"];
Print[diff];
Print[""];
If[Simplify[diff] === 0,
  Print["✓ GENERAL EXPRESSION MATCHES"],
  Print["✗ MISMATCH! Simplifying difference:"];
  Print[Simplify[diff]]
];

(* ============================================================ *)
Print[""];
Print["===== PART 2: S-channel specialization (3,3,4,4) ====="];
Print[""];

(* Set vi=vj=v3, vk=vl=v4, with v3^2=v4^2=0, v3.v4=v34 *)
ScalarProduct[v3, v3] = 0;
ScalarProduct[v4, v4] = 0;

(* FeynCalc direct S-channel *)
VfullS = FV[v3, LorentzIndex[nu1]] * FV[v3, LorentzIndex[nu2]] *
         FV[v4, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VcontractedS = ExpandScalarProduct[Contract[VfullS]];
Print["FeynCalc S-channel (v3,v3,v4,v4):"];
Print[VcontractedS];
Print[""];

(* Document's S-channel: substitute vi=vj=v3, vk=vl=v4 *)
VdocSubS = VdocExpanded /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4};
(* Expand and apply v3^2=v4^2=0 *)
VdocSubS = ExpandScalarProduct[VdocSubS];
Print["Document S-channel (substitute vi=vj=v3, vk=vl=v4):"];
Print[VdocSubS];
Print[""];

diffS = ExpandScalarProduct[VcontractedS - VdocSubS];
Print["Difference for S-channel:"];
Print[diffS];
Print[""];
If[Simplify[diffS] === 0,
  Print["✓ S-CHANNEL MATCHES"],
  Print["✗ S-CHANNEL MISMATCH! Simplifying:"];
  Print[Simplify[diffS]]
];

(* Expected S-channel result from derivation doc: z (v3.(l1-l2)) (v4.(l3-l4)) *)
(* Note: v3.v4 = v34 = z *)
expectedS = Pair[v3, v4] * (Pair[v3, l1] - Pair[v3, l2]) * (Pair[v4, l3] - Pair[v4, l4]);
expectedS = ExpandScalarProduct[expectedS];
Print[""];
Print["Expected S-channel: z (v3.(l1-l2)) (v4.(l3-l4)) ="];
Print[expectedS];
diffSexp = ExpandScalarProduct[VcontractedS - expectedS];
Print[""];
Print["FeynCalc S-channel vs expected:"];
If[Simplify[diffSexp] === 0,
  Print["✓ S-CHANNEL MATCHES EXPECTED"],
  Print["✗ S-CHANNEL vs EXPECTED MISMATCH:"];
  Print[Simplify[diffSexp]]
];

(* ============================================================ *)
Print[""];
Print["===== PART 3: T-channel specialization (3,4,3,4) ====="];
Print[""];

(* FeynCalc direct T-channel: vi=v3, vj=v4, vk=v3, vl=v4 *)
VfullT = FV[v3, LorentzIndex[nu1]] * FV[v4, LorentzIndex[nu2]] *
         FV[v3, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VcontractedT = ExpandScalarProduct[Contract[VfullT]];
Print["FeynCalc T-channel (v3,v4,v3,v4):"];
Print[VcontractedT];
Print[""];

(* Document's T-channel specialization: vi=v3, vj=v4, vk=v3, vl=v4 *)
VdocSubT = VdocExpanded /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4};
VdocSubT = ExpandScalarProduct[VdocSubT];
Print["Document T-channel (substitute vi=v3, vj=v4, vk=v3, vl=v4):"];
Print[VdocSubT];
Print[""];

diffT = ExpandScalarProduct[VcontractedT - VdocSubT];
Print["Difference for T-channel:"];
If[Simplify[diffT] === 0,
  Print["✓ T-CHANNEL MATCHES"],
  Print["✗ T-CHANNEL MISMATCH! Simplifying:"];
  Print[Simplify[diffT]]
];

(* ============================================================ *)
Print[""];
Print["===== PART 4: Key assertions check ====="];
Print[""];

(* Check: S-channel has only Type III terms (G1, G2, G3 should vanish) *)
Print["Checking S-channel degeneration:"];
Print["G1 (should vanish: vi.vj=0): ", Simplify[G1 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}]];
Print["G2 (should vanish: vi.vj=0): ", Simplify[G2 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}]];
Print["G3 (should vanish: vk.vl=0): ", Simplify[G3 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}]];
Print["G4a (survives): ", ExpandScalarProduct[G4a /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}]];
Print["G4b (survives): ", ExpandScalarProduct[G4b /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}]];
Print["G4c (survives): ", ExpandScalarProduct[G4c /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}]];
Print["G4d (survives): ", ExpandScalarProduct[G4d /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}]];

Print[""];
Print["===== VERIFICATION COMPLETE ====="];