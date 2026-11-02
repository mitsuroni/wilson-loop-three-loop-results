(* verify_pnum_structure.wl -- Explore Pnum structure for 4e-3g+3g *)

$FeynCalcStartupMessages = False;
$FCCheckContext = False;
Needs["FeynCalc`"];

sp[a_, b_] := Pair[Momentum[a], Momentum[b]];

Print["===== PNUM STRUCTURE ANALYSIS ====="];
Print[""];

(* ============================================================ *)
(* Part 1: S-channel Pnum structure *)
(* ============================================================ *)

Print["Part 1: S-channel Pnum decomposition"];

(* O_S = -z * (v3·∂1 - v3·∂2) * (v4·∂3 - v4·∂4) *)
(* = -z * [v3·∂1 v4·∂3 - v3·∂1 v4·∂4 - v3·∂2 v4·∂3 + v3·∂2 v4·∂4] *)

(* Each term (v3·∂i)(v4·∂j) acts on B_5^n: *)
(* (v3·∂i)(v4·∂j) B_5^n = n(n-1) B_5^{n-2} F_i G_j + n B_5^{n-1} H_{ij} *)

(* where F_i = v3·∂i B_5, G_j = v4·∂j B_5, H_{ij} = v3·∂i G_j *)

(* For S-channel with v3·X1 = v3·X2 = -z, v4·X3 = v4·X4 = 0: *)

(* F_1 = v3·∂1 B_5 = -2*a1*[a2*SB*0 + a5*a3*(v3·X1 - v3·X3) + a5*a4*(v3·X1 - v3·X4)] *)
(* v3·X1 = -z, v3·X3 = z*(s41 - 1), v3·X4 = z*(s42 - 1) *)
(* F_1 = -2*a1*[a5*a3*(-z - z*(s41-1)) + a5*a4*(-z - z*(s42-1))] *)
(*     = -2*a1*[a5*a3*(-z*s41) + a5*a4*(-z*s42)] *)
(*     = 2*a1*a5*z*(a3*s41 + a4*s42) *)

(* F_2 = v3·∂2 B_5 = -2*a2*[a1*SB*0 + a5*a3*(v3·X2 - v3·X3) + a5*a4*(v3·X2 - v3·X4)] *)
(*     = 2*a2*a5*z*(a3*s41 + a4*s42) *)

(* So F_1 - F_2 = 2*a5*z*(a1-a2)*(a3*s41 + a4*s42) = F (the doc's F) *)

(* G_3 = v4·∂3 B_5 = -2*a3*[a4*SA*0 + a5*a1*(v4·X3 - v4·X1) + a5*a2*(v4·X3 - v4·X2)] *)
(* v4·X3 = 0, v4·X1 = z*(s31-1), v4·X2 = z*(s32-1) *)
(* G_3 = -2*a3*[a5*a1*(0 - z*(s31-1)) + a5*a2*(0 - z*(s32-1))] *)
(*     = -2*a3*[a5*a1*z*(1-s31) + a5*a2*z*(1-s32)] *)
(*     = -2*a3*a5*z*(a1*(1-s31) + a2*(1-s32)) *)

(* G_4 = v4·∂4 B_5 = -2*a4*[a3*SA*0 + a5*a1*(v4·X4 - v4·X1) + a5*a2*(v4·X4 - v4·X2)] *)
(* v4·X4 = 0 *)
(* G_4 = -2*a4*a5*z*(a1*(1-s31) + a2*(1-s32)) *)

(* So G_3 - G_4 = -2*a5*z*(a3-a4)*(a1*(1-s31) + a2*(1-s32)) = G (the doc's G) *)

(* H_{ij} = v3·∂i G_j *)
(* For i=1, j=3: *)
(* G_3 depends on X1 through v4·X1. v3·∂1(v4·X1) = v3·v4 = z. *)
(* H_{13} = v3·∂1[-2*a3*a5*z*a1*(1-s31)] = -2*a3*a5*z*a1 * v3·∂1(1-s31) *)
(* 1-s31 = -v4·X1/z. v3·∂1(-v4·X1/z) = -v3·v4/z = -z*v_ref/z = -v_ref? *)

(* Hmm, let me be more careful. s31 = x1*x2 in hypercube. X1 = v1+v2 + s31*v3. *)
(* v4·X1 = v4·(v1+v2) + s31*v4·v3 = v4·(v1+v2) + s31*z*v_ref *)
(* Wait, v4·v1 = v14 = v_ref, v4·v2 = v24 = -(1+z)*v_ref. *)
(* v4·(v1+v2) = v_ref - (1+z)*v_ref = -z*v_ref. *)
(* So v4·X1 = -z*v_ref + s31*z*v_ref = z*v_ref*(s31-1). *)
(* 1-s31 = -v4·X1/(z*v_ref). *)

(* v3·∂1(v4·X1) = v3·v4 = z*v_ref. *)
(* v3·∂1(1-s31) = v3·∂1(-v4·X1/(z*v_ref)) = -v3·v4/(z*v_ref) = -z*v_ref/(z*v_ref) = -1. *)

(* Wait, that can't be right. ∂1 is the derivative with respect to X1, not s31. *)
(* v3·∂1 is the directional derivative along v3. *)
(* v3·∂1(v4·X1) = v3·(v4) = v3·v4 = z*v_ref. ✓ *)
(* But how does v3·∂1 act on s31? s31 is a parameter, not a function of X1. *)
(* Actually, X1 = V1 + s31*v3, so s31 = (X1 - V1)·v3 / v3^2... but v3^2 = 0! *)

(* This is the key subtlety. On a lightlike edge, s31 is not uniquely determined by X1. *)
(* The derivative v3·∂1 acts on functions of X1, and s31 is a parameter that labels X1. *)
(* We need to use the chain rule: v3·∂1 f(X1(s31)) = v3·(∂X1/∂s31) ∂f/∂s31 = v3·v3 ∂f/∂s31 = 0! *)

(* So v3·∂1 f(s31) = 0 for any function of s31 alone! *)
(* But G_3 depends on X1 through v4·X1, not through s31. *)
(* v4·X1 is a genuine D-dimensional function of X1. *)
(* v3·∂1(v4·X1) = v3·v4 = z*v_ref ≠ 0. *)

(* Let me redo this more carefully. *)
(* G_3 = -2*a3*a5*z*(a1*(1-s31) + a2*(1-s32)) *)
(* Here s31 and s32 are the Schwinger parameters. *)
(* 1-s31 = -v4·X1/(z*v_ref). This is a function of X1. *)

(* v3·∂1 G_3 = v3·∂1[-2*a3*a5*z*a1*(-v4·X1/(z*v_ref))] *)
(*            = v3·∂1[2*a3*a5*a1*v4·X1/v_ref] *)
(*            = 2*a3*a5*a1*(v3·v4)/v_ref *)
(*            = 2*a3*a5*a1*z*v_ref/v_ref *)
(*            = 2*a3*a5*a1*z *)

(* So H_{13} = 2*a3*a5*a1*z *)

(* H_{14} = v3·∂1 G_4 = v3·∂1[-2*a4*a5*z*a1*(1-s31)] = 2*a4*a5*a1*z *)
(* H_{23} = v3·∂2 G_3 = 2*a3*a5*a2*z *)
(* H_{24} = v3·∂2 G_4 = 2*a4*a5*a2*z *)

(* H_13 - H_14 - H_23 + H_24 = 2*a5*z*[a3*a1 - a4*a1 - a3*a2 + a4*a2] *)
(*                            = 2*a5*z*(a1-a2)*(a3-a4) *)
(* This matches the doc's H = 2*a5*z*(a3-a4)*(a1-a2)! ✓ *)

(* Now, let's assemble the Pnum for the S-channel. *)
(* O_S = -z * [T13 - T14 - T23 + T24] where Tij = (v3·∂i)(v4·∂j) *)
(* Each Tij contributes: n(n-1)B_5^{n-2} F_i G_j + n B_5^{n-1} H_{ij} *)

(* Pnum_FG = -z * [F_1*G_3 - F_1*G_4 - F_2*G_3 + F_2*G_4] *)
(*         = -z * (F_1 - F_2)*(G_3 - G_4) *)
(*         = -z * F * G *)

(* Pnum_H = -z * [H_13 - H_14 - H_23 + H_24] *)
(*        = -z * H *)

(* O_S(B_5^n) = n(n-1) B_5^{n-2} * (-z*F*G) + n B_5^{n-1} * (-z*H) *)
(*            = -z*[n(n-1) B_5^{n-2} F*G + n B_5^{n-1} H] *)

Print["S-channel Pnum structure:"];
Print["Pnum = -z * [F*G * (factor) + H * (factor)]"];
Print["where F = v3·(∂1-∂2)B_5, G = v4·(∂3-∂4)B_5, H = v3·(∂1-∂2)G"];
Print[""];
Print["F = 2*a5*z*(a1-a2)*(a3*s41 + a4*s42)"];
Print["G = 2*a5*z*(a4-a3)*(a1*(1-s31) + a2*(1-s32))"];
Print["H = 2*a5*z*(a3-a4)*(a1-a2)"];
Print[""];
Print["Note: F*G ∝ z^2, H ∝ z. Both ∝ a5 or a5^2."];
Print["The prefactor -z makes the overall Pnum ∝ z^3 or z^2."];

(* ============================================================ *)
(* Part 2: Can we write Pnum like (1/4) sum α_i α_j W_ij? *)
(* ============================================================ *)

Print[""];
Print["Part 2: Comparison with TGV structure"];

Print["TGV (3e-3g):  Pnum = (1/4)(α1α2 W12 + α2α3 W23 + α3α1 W31)"];
Print["  - Single B_5 power: S^{-(2-2ε)}"];
Print["  - O is linear → single term from chain rule"];
Print["  - W_ij are simple polynomials in x_i, z"];
Print[""];
Print["4e-3g+3g (S-channel):"];
Print["  Pnum = α5^2 * Pnum0 * (2-3ε) / B_5^{3-3ε}  +  α5*(α3-α4)*(α1-α2) / B_5^{2-3ε}"];
Print["  - TWO B_5 powers: -(3-3ε) and -(2-3ε)"];
Print["  - O is quadratic → two terms from chain rule (FG + H)"];
Print["  - Extra U^{-(1+2ε)} in denominator (from 2-loop topology)"];

Print[""];
Print["Key structural differences:"];
Print["1. TWO terms with different B_5 powers (unavoidable — quadratic O)"];
Print["2. U polynomial in denominator (unavoidable — 2-loop topology)"];
Print["3. Pnum depends on α5 (the bridge Feynman parameter) prominently"];
Print["4. No simple (1/4)∑α_iα_j W_{ij} factorization — the FG and H terms"];
Print["   have different α-structure and different B_5 scaling"];

(* ============================================================ *)
(* Part 3: General Pnum structure *)
(* ============================================================ *)

Print[""];
Print["Part 3: General Pnum structure for arbitrary (a,b,c,d)"];

Print["O = sum over 15 terms, each of form (v·v prefactor) × (v_p·∂_i)(v_q·∂_j)"];
Print[""];
Print["For each term: (v_p·∂_i)(v_q·∂_j) B_5^n"];
Print["  = n(n-1) B_5^{n-2} F_i G_j + n B_5^{n-1} H_{ij}"];
Print[""];
Print["where:"];
Print["  F_i = v_p·∂_i B_5 = -2α_i ∑_{j≠i} α_j coeff_{ij} v_p·(X_i - X_j)"];
Print["  G_j = v_q·∂_j B_5 = -2α_j ∑_{k≠j} α_k coeff_{jk} v_q·(X_j - X_k)"];
Print["  H_{ij} = v_p·∂_i G_j = -2α_i α_j coeff_{ij} (v_p·v_q) × (sign)"];
Print[""];
Print["coeff_{ij} = Σ_B (i,j∈{1,2}), Σ_A (i,j∈{3,4}), α_5 (cross-group)"];
Print[""];
Print["General Pnum = Σ_{O terms} (v·v prefactor) × [n(n-1) F_i G_j + n H_{ij}]"];
Print[""];
Print["n = -(1-3ε), n(n-1) = (1-3ε)(2-3ε)"];

(* ============================================================ *)
Print[""];
Print["===== ANALYSIS COMPLETE ====="];