(* verify_v2.wl -- Fixed verification with consistent FeynCalc notation *)

$FeynCalcStartupMessages = False;
$FCCheckContext = False;
Needs["FeynCalc`"];

Print["===== VERIFICATION OF 4e-3g+3g GENERAL DERIVATION (v2) ====="];
Print[""];

(* ============================================================ *)
(* PART 1: General V_{2TGV} - compare FeynCalc vs document *)
(* ============================================================ *)

(* TGV kinematic structures *)
GvA[nu1_, nu2_, nu5_] :=
  MT[nu1, nu2] * FV[l2 - l1, nu5] +
  MT[nu2, nu5] * FV[-l1 - 2*l2, nu1] +
  MT[nu5, nu1] * FV[2*l1 + l2, nu2];

GvB[nu3_, nu4_, nu6_] :=
  MT[nu3, nu4] * FV[l4 - l3, nu6] +
  MT[nu4, nu6] * FV[-l1 - l2 - l4, nu3] +
  MT[nu6, nu3] * FV[l1 + l2 + l3, nu4];

(* FeynCalc general contraction *)
VfullFC = FV[vi, LorentzIndex[nu1]] * FV[vj, LorentzIndex[nu2]] *
          FV[vk, LorentzIndex[nu3]] * FV[vl, LorentzIndex[nu4]] *
          MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
          GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
          GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

Vfc = Contract[VfullFC] // ExpandScalarProduct;

(* Document's expression using FeynCalc SP (ScalarProduct) *)
(* G1: (vi.vj)(vk.vl) (l2-l1).(l4-l3) *)
G1 = SP[vi, vj] * SP[vk, vl] *
     (SP[l2, l4] - SP[l2, l3] - SP[l1, l4] + SP[l1, l3]);

(* G2: (vi.vj)[...] *)
G2 = SP[vi, vj] * (
  2*SP[vk, l1]*SP[vl, l2] - 2*SP[vk, l2]*SP[vl, l1] +
  (SP[vk, l2] - SP[vk, l1])*SP[vl, l3] +
  (SP[vl, l1] - SP[vl, l2])*SP[vk, l4]
);

(* G3: 2(vk.vl)[...] *)
G3 = 2*SP[vk, vl] * (
  (-SP[vi, l1] - 2*SP[vi, l2])*(SP[vj, l4] - SP[vj, l3]) +
  (2*SP[vj, l1] + SP[vj, l2])*(SP[vi, l4] - SP[vi, l3])
);

(* G4a-d *)
G4a = SP[vi, vk] * (2*SP[vj, l1] + SP[vj, l2]) *
      (SP[vl, l3] - SP[vl, l1] - SP[vl, l2]);
G4b = SP[vi, vl] * (2*SP[vj, l1] + SP[vj, l2]) *
      (SP[vk, l1] + SP[vk, l2] - SP[vk, l4]);
G4c = -SP[vj, vk] * (SP[vi, l1] + 2*SP[vi, l2]) *
      (SP[vl, l3] - SP[vl, l1] - SP[vl, l2]);
G4d = -SP[vj, vl] * (SP[vi, l1] + 2*SP[vi, l2]) *
      (SP[vk, l1] + SP[vk, l2] - SP[vk, l4]);

Vdoc = (G1 + G2 + G3 + G4a + G4b + G4c + G4d) // ExpandScalarProduct;

Print["Part 1: General expression comparison"];
diff1 = Vfc - Vdoc // Simplify;
If[diff1 === 0,
  Print["✓ GENERAL EXPRESSION MATCHES"],
  Print["✗ MISMATCH. Difference:"];
  Print[diff1]
];

(* ============================================================ *)
(* PART 2: S-channel (3,3,4,4) *)
(* ============================================================ *)

ScalarProduct[v3, v3] = 0;
ScalarProduct[v4, v4] = 0;

(* FeynCalc direct *)
VfullS = FV[v3, LorentzIndex[nu1]] * FV[v3, LorentzIndex[nu2]] *
         FV[v4, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcS = Contract[VfullS] // ExpandScalarProduct;

(* Document's S-channel: substitute vi=vj=v3, vk=vl=v4 *)
VdocS = (Vdoc /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;

(* Expected S-channel: z (v3.(l1-l2)) (v4.(l3-l4)) *)
expectedS = SP[v3, v4] * (SP[v3, l1] - SP[v3, l2]) * (SP[v4, l3] - SP[v4, l4]);

Print[""];
Print["Part 2: S-channel specialization (3,3,4,4)"];

diffS_doc = VfcS - VdocS // Simplify;
If[diffS_doc === 0,
  Print["✓ S-CHANNEL: Document substitution matches FeynCalc"],
  Print["✗ S-CHANNEL: Document substitution MISMATCH:"];
  Print[diffS_doc]
];

diffS_exp = VfcS - expectedS // Simplify;
If[diffS_exp === 0,
  Print["✓ S-CHANNEL: FeynCalc matches expected z(v3.(l1-l2))(v4.(l3-l4))"],
  Print["✗ S-CHANNEL: FeynCalc vs expected MISMATCH:"];
  Print[diffS_exp]
];

(* ============================================================ *)
(* PART 3: Check S-channel degeneration *)
(* ============================================================ *)

Print[""];
Print["Part 3: S-channel degeneration check"];

(* Each group after substitution *)
G1s = (G1 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;
G2s = (G2 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;
G3s = (G3 /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;
G4as = (G4a /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;
G4bs = (G4b /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;
G4cs = (G4c /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;
G4ds = (G4d /. {vi -> v3, vj -> v3, vk -> v4, vl -> v4}) // ExpandScalarProduct;

checks = {
  {"G1 (v3.v3=0 → vanish)", G1s, 0},
  {"G2 (v3.v3=0 → vanish)", G2s, 0},
  {"G3 (v4.v4=0 → vanish)", G3s, 0},
  {"G4a+G4b+G4c+G4d (survive)", G4as+G4bs+G4cs+G4ds, expectedS}
};

Do[
  diff = check[[2]] - check[[3]] // Simplify;
  If[diff === 0,
    Print["✓ ", check[[1]]],
    Print["✗ ", check[[1]], " — diff: ", diff]
  ],
  {check, checks}
];

(* ============================================================ *)
(* PART 4: T-channel (3,4,3,4) *)
(* ============================================================ *)

VfullT = FV[v3, LorentzIndex[nu1]] * FV[v4, LorentzIndex[nu2]] *
         FV[v3, LorentzIndex[nu3]] * FV[v4, LorentzIndex[nu4]] *
         MT[LorentzIndex[nu5], LorentzIndex[nu6]] *
         GvA[LorentzIndex[nu1], LorentzIndex[nu2], LorentzIndex[nu5]] *
         GvB[LorentzIndex[nu3], LorentzIndex[nu4], LorentzIndex[nu6]];

VfcT = Contract[VfullT] // ExpandScalarProduct;

VdocT = (Vdoc /. {vi -> v3, vj -> v4, vk -> v3, vl -> v4}) // ExpandScalarProduct;

Print[""];
Print["Part 4: T-channel specialization (3,4,3,4)"];

diffT = VfcT - VdocT // Simplify;
If[diffT === 0,
  Print["✓ T-CHANNEL: Document substitution matches FeynCalc"],
  Print["✗ T-CHANNEL: Document substitution MISMATCH:"];
  Print[diffT]
];

(* ============================================================ *)
Print[""];
Print["===== VERIFICATION COMPLETE ====="];