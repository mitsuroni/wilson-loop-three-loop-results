(* verify_O_operator.wl -- Verify O(a,b,c,d) = L with l->i*dX *)

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
(* PART 1: Generate O from FeynCalc V_{2TGV} by l->i*dX *)
(* ============================================================ *)

Print["===== PART 1: O from FeynCalc (l->i*dX) ====="];

(* Compute FeynCalc V_{2TGV} *)
VfullFC = FV[vi, LorentzIndex[nu1]] * FV[vj, LorentzIndex[nu2]] *
          FV[vk, LorentzIndex[nu3]] * FV[vl, LorentzIndex[nu4]] *
          MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
          GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
          GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

Vfc = Contract[VfullFC] // ExpandScalarProduct;

(* Replace l_i -> i * dX_i *)
(* sp[l_i, v_j] -> i * sp[dX_i, v_j] *)
(* sp[l_i, l_j] -> i^2 * sp[dX_i, dX_j] = -sp[dX_i, dX_j] *)

OfromFC = Vfc /. {
  sp[l1, v_] :> I * sp[dX1, v],
  sp[l2, v_] :> I * sp[dX2, v],
  sp[l3, v_] :> I * sp[dX3, v],
  sp[l4, v_] :> I * sp[dX4, v],
  sp[v_, l1] :> I * sp[v, dX1],
  sp[v_, l2] :> I * sp[v, dX2],
  sp[v_, l3] :> I * sp[v, dX3],
  sp[v_, l4] :> I * sp[v, dX4],
  sp[l1, l1] :> -sp[dX1, dX1],
  sp[l2, l2] :> -sp[dX2, dX2],
  sp[l3, l3] :> -sp[dX3, dX3],
  sp[l4, l4] :> -sp[dX4, dX4],
  sp[l1, l2] :> -sp[dX1, dX2],
  sp[l1, l3] :> -sp[dX1, dX3],
  sp[l1, l4] :> -sp[dX1, dX4],
  sp[l2, l1] :> -sp[dX2, dX1],
  sp[l2, l3] :> -sp[dX2, dX3],
  sp[l2, l4] :> -sp[dX2, dX4],
  sp[l3, l1] :> -sp[dX3, dX1],
  sp[l3, l2] :> -sp[dX3, dX2],
  sp[l3, l4] :> -sp[dX3, dX4],
  sp[l4, l1] :> -sp[dX4, dX1],
  sp[l4, l2] :> -sp[dX4, dX2],
  sp[l4, l3] :> -sp[dX4, dX3]
};

OfromFC = OfromFC // ExpandScalarProduct;

(* ============================================================ *)
(* PART 2: Document's O operator *)
(* ============================================================ *)

(* O_I: -(vi.vj)(vk.vl) (dX2-dX1).(dX4-dX3) *)
OI = -sp[vi, vj] * sp[vk, vl] *
     (sp[dX2, dX4] - sp[dX2, dX3] - sp[dX1, dX4] + sp[dX1, dX3]);

(* O_II^{(ab)-G2}: -v_ab [...] *)
OIIabG2 = -sp[vi, vj] * (
  2*sp[vc, dX1]*sp[vd, dX2] - 2*sp[vc, dX2]*sp[vd, dX1] +
  (sp[vc, dX2] - sp[vc, dX1])*sp[vd, dX3] +
  (sp[vd, dX1] - sp[vd, dX2])*sp[vc, dX4]
);

(* O_II^{(ab)-G3}: -4*v_ab [(vd.dX1)(vc.dX2) - (vc.dX1)(vd.dX2)] *)
OIIabG3 = -4*sp[vi, vj] * (sp[vd, dX1]*sp[vc, dX2] - sp[vc, dX1]*sp[vd, dX2]);

(* O_II^{(cd)}: -v_cd [(va.dX1+2va.dX2)(vb.dX3-vb.dX4) - (2vb.dX1+vb.dX2)(va.dX3-va.dX4)] *)
OIIcd = -sp[vk, vl] * (
  (sp[vi, dX1] + 2*sp[vi, dX2])*(sp[vj, dX3] - sp[vj, dX4]) -
  (2*sp[vj, dX1] + sp[vj, dX2])*(sp[vi, dX3] - sp[vi, dX4])
);

(* O_III^{(G3)}: 4 cross terms *)
OIIIg3a = -2*sp[vi, vk] * (2*sp[vj, dX1] + sp[vj, dX2]) * (sp[vd, dX1] + sp[vd, dX2]);
OIIIg3b = 2*sp[vi, vl] * (2*sp[vj, dX1] + sp[vj, dX2]) * (sp[vc, dX1] + sp[vc, dX2]);
OIIIg3c = 2*sp[vj, vk] * (sp[vi, dX1] + 2*sp[vi, dX2]) * (sp[vd, dX1] + sp[vd, dX2]);
OIIIg3d = -2*sp[vj, vl] * (sp[vi, dX1] + 2*sp[vi, dX2]) * (sp[vc, dX1] + sp[vc, dX2]);
OIIIg3 = OIIIg3a + OIIIg3b + OIIIg3c + OIIIg3d;

(* O_III^{(G4)}: 4 original terms *)
OIIIg4a = -sp[vi, vk] * (2*sp[vj, dX1] + sp[vj, dX2]) * (sp[vd, dX3] - sp[vd, dX1] - sp[vd, dX2]);
OIIIg4b = -sp[vi, vl] * (2*sp[vj, dX1] + sp[vj, dX2]) * (sp[vc, dX1] + sp[vc, dX2] - sp[vc, dX4]);
OIIIg4c = sp[vj, vk] * (sp[vi, dX1] + 2*sp[vi, dX2]) * (sp[vd, dX3] - sp[vd, dX1] - sp[vd, dX2]);
OIIIg4d = sp[vj, vl] * (sp[vi, dX1] + 2*sp[vi, dX2]) * (sp[vc, dX1] + sp[vc, dX2] - sp[vc, dX4]);
OIIIg4 = OIIIg4a + OIIIg4b + OIIIg4c + OIIIg4d;

Odoc = (OI + OIIabG2 + OIIabG3 + OIIcd + OIIIg3 + OIIIg4) // ExpandScalarProduct;

diff = OfromFC - Odoc // Simplify;
Print["O(FeynCalc) - O(Document) = "];
If[diff === 0,
  Print["✓ ZERO — O operator matches"],
  Print["✗ NON-ZERO. Difference:"];
  Print[diff]
];

(* ============================================================ *)
(* PART 3: S-channel O reduction *)
(* ============================================================ *)

Print[""];
Print["===== PART 2: S-channel O reduction ====="];

ScalarProduct[v3, v3] = 0;
ScalarProduct[v4, v4] = 0;

(* FeynCalc S-channel l->i*dX *)
VfullS = FV[v3, LorentzIndex[nu1]] * FV[v3, LorentzIndex[nu2]] *
         FV[v4, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcS = Contract[VfullS] // ExpandScalarProduct;

OfromFCS = VfcS /. {
  sp[l1, v_] :> I * sp[dX1, v],
  sp[l2, v_] :> I * sp[dX2, v],
  sp[l3, v_] :> I * sp[dX3, v],
  sp[l4, v_] :> I * sp[dX4, v],
  sp[v_, l1] :> I * sp[v, dX1],
  sp[v_, l2] :> I * sp[v, dX2],
  sp[v_, l3] :> I * sp[v, dX3],
  sp[v_, l4] :> I * sp[v, dX4],
  sp[l1, l2] :> -sp[dX1, dX2],
  sp[l1, l3] :> -sp[dX1, dX3],
  sp[l1, l4] :> -sp[dX1, dX4],
  sp[l2, l3] :> -sp[dX2, dX3],
  sp[l2, l4] :> -sp[dX2, dX4],
  sp[l3, l4] :> -sp[dX3, dX4]
};

OfromFCS = OfromFCS // ExpandScalarProduct // Simplify;

(* Document's O specialized to S-channel *)
OdocS = (Odoc /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4, vc -> v4, vd -> v4}) // ExpandScalarProduct // Simplify;

(* Expected O_S = -z (v3.(dX1-dX2)) (v4.(dX3-dX4)) *)
(* Wait, let me check: v3.v4 = v34 = z *)
expectedOS = -sp[v3, v4] * (sp[v3, dX1] - sp[v3, dX2]) * (sp[v4, dX3] - sp[v4, dX4]);

Print["O_S from FeynCalc: ", OfromFCS];
Print[""];
Print["Expected O_S = -z(v3.(dX1-dX2))(v4.(dX3-dX4)): ", expectedOS];

diffS_fc = OfromFCS - expectedOS // Simplify;
diffS_doc = OdocS - expectedOS // Simplify;

Print[""];
Print["FeynCalc O_S vs expected: ", If[diffS_fc===0, "✓", "✗ " <> ToString[diffS_fc]]];
Print["Document O_S vs expected: ", If[diffS_doc===0, "✓", "✗ " <> ToString[diffS_doc]]];

(* ============================================================ *)
(* PART 4: T-channel O reduction *)
(* ============================================================ *)

Print[""];
Print["===== PART 3: T-channel O reduction ====="];

VfullT = FV[v3, LorentzIndex[nu1]] * FV[v4, LorentzIndex[nu2]] *
         FV[v3, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcT = Contract[VfullT] // ExpandScalarProduct;

OfromFCT = VfcT /. {
  sp[l1, v_] :> I * sp[dX1, v],
  sp[l2, v_] :> I * sp[dX2, v],
  sp[l3, v_] :> I * sp[dX3, v],
  sp[l4, v_] :> I * sp[dX4, v],
  sp[v_, l1] :> I * sp[v, dX1],
  sp[v_, l2] :> I * sp[v, dX2],
  sp[v_, l3] :> I * sp[v, dX3],
  sp[v_, l4] :> I * sp[v, dX4],
  sp[l1, l2] :> -sp[dX1, dX2],
  sp[l1, l3] :> -sp[dX1, dX3],
  sp[l1, l4] :> -sp[dX1, dX4],
  sp[l2, l3] :> -sp[dX2, dX3],
  sp[l2, l4] :> -sp[dX2, dX4],
  sp[l3, l4] :> -sp[dX3, dX4]
};

OfromFCT = OfromFCT // ExpandScalarProduct // Simplify;

OdocT = (Odoc /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4, vc -> v3, vd -> v4}) // ExpandScalarProduct // Simplify;

diffT = OfromFCT - OdocT // Simplify;
Print["O_T from FeynCalc vs Document: ", If[diffT===0, "✓ MATCHES", "✗ " <> ToString[diffT]]];

Print[""];
Print["===== VERIFICATION COMPLETE ====="];