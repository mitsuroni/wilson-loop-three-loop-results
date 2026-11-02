(* verify_symanzik_5pt.wl -- Step-by-step verification of 5-point Symanzik parametrization *)

$FeynCalcStartupMessages = False;
$FCCheckContext = False;
Needs["FeynCalc`"];

Print["===== 5-POINT SYMANZIK PARAMETRIZATION VERIFICATION ====="];
Print[""];

(* ============================================================ *)
(* Step 1: Define the 5 propagators and Q *)
(* ============================================================ *)

Print["Step 1: Define Q = sum(alpha_i * Delta_i)"];

(* Delta_i = -(X_i - z)^2 + i0. We drop i0 for algebraic verification. *)
(* Delta_1 = -(X1 - zA)^2, Delta_2 = -(X2 - zA)^2, Delta_5 = -(zA - zB)^2 *)
(* Delta_3 = -(X3 - zB)^2, Delta_4 = -(X4 - zB)^2 *)

(* Use FeynCalc scalar products for convenience *)
sp[a_, b_] := Pair[Momentum[a], Momentum[b]];

(* For now, treat X_i and zA, zB as vectors *)
(* We'll work with symbolic squared distances *)

(* Define Delta_ij explicitly *)
(* Delta_1 = -(X1 - zA)^2 = -X1^2 + 2 X1·zA - zA^2 *)
(* Similarly for others *)

(* Let's define Q in terms of scalar products *)
(* Q = a1*( -X1^2 + 2*X1·zA - zA^2 ) + a2*( -X2^2 + 2*X2·zA - zA^2 )
      + a5*( -zA^2 + 2*zA·zB - zB^2 )
      + a3*( -X3^2 + 2*X3·zB - zB^2 ) + a4*( -X4^2 + 2*X4·zB - zB^2 ) *)

(* ============================================================ *)
(* Step 2: Extract Sigma, J, C from Q *)
(* ============================================================ *)

Print["Step 2: Extract Sigma matrix, J vector, C scalar"];

(* Q = z^T (Sigma ⊗ I_D) z - 2 J·z + C *)
(* where z = (zA, zB)^T *)

(* zA^2 coefficient: -(a1+a2+a5) → Sigma_11 = a1+a2+a5 *)
(* zB^2 coefficient: -(a3+a4+a5) → Sigma_22 = a3+a4+a5 *)
(* zA·zB coefficient: 2*a5 → Sigma_12 = Sigma_21 = -a5 (with the -2*J·z convention) *)

(* Let's verify: *)
(* Q = -(a1+a2+a5)*zA^2 - (a3+a4+a5)*zB^2 + 2*a5*zA·zB
      + 2*a1*X1·zA + 2*a2*X2·zA - 2*a3*X3·zB - 2*a4*X4·zB
      - a1*X1^2 - a2*X2^2 - a3*X3^2 - a4*X4^2 *)

(* Wait, there's a sign issue. Let me be more careful. *)
(* Delta_i = -(X_i - z)^2 = -(X_i^2 - 2X_i·z + z^2) = -X_i^2 + 2X_i·z - z^2 *)
(* For Delta_5 = -(zA - zB)^2 = -zA^2 + 2zA·zB - zB^2 *)

(* So Q = a1*(-X1^2 + 2X1·zA - zA^2) + a2*(-X2^2 + 2X2·zA - zA^2)
         + a5*(-zA^2 + 2zA·zB - zB^2)
         + a3*(-X3^2 + 2X3·zB - zB^2) + a4*(-X4^2 + 2X4·zB - zB^2) *)

(* Quadratic terms in z: *)
(* zA^2: -(a1+a2+a5) = -Sigma_11  → Sigma_11 = a1+a2+a5 *)
(* zB^2: -(a3+a4+a5) = -Sigma_22  → Sigma_22 = a3+a4+a5 *)
(* zA·zB: 2*a5 → this comes from -2*Sigma_12*zA·zB in the form z^T Sigma z *)
(* For z^T Sigma z = Sigma_11*zA^2 + Sigma_22*zB^2 + 2*Sigma_12*zA·zB *)
(* We have -Sigma_11*zA^2 - Sigma_22*zB^2 + 2*a5*zA·zB *)
(* But Q = z^T Sigma z - 2J·z + C, and z^T Sigma z = Sigma_11*zA^2 + 2*Sigma_12*zA·zB + Sigma_22*zB^2 *)
(* So Q = Sigma_11*zA^2 + 2*Sigma_12*zA·zB + Sigma_22*zB^2 - 2J·z + C *)
(* We have Q = -(a1+a2+a5)*zA^2 + 2*a5*zA·zB - (a3+a4+a5)*zB^2 + linear + constant *)
(* So: Sigma_11 = -(a1+a2+a5)? No, that doesn't work. *)

(* Hmm, let me reconsider. The standard form is Q = -[z^T (Sigma ⊗ I_D) z - 2J·z + C] *)
(* OR Q = z^T (Sigma ⊗ I_D) z - 2J·z + C *)
(* The sign convention depends on whether Delta includes the minus sign. *)

(* Let me just compute Q explicitly and then extract the coefficients. *)

(* The Q form is: Q = alpha_i * Delta_i *)
(* After Feynman parametrization, we have Q^{-5a} in the denominator *)
(* The Gaussian integral formula is: ∫ d^{2D}z Q^{-N} = (standard formula) *)

(* Let me just verify the matrix form by direct computation. *)

(* ============================================================ *)
(* Step 3: Verify Sigma and adj(Sigma) *)
(* ============================================================ *)

Print["Step 3: Verify Sigma, det(Sigma), adj(Sigma)"];

(* Sigma = {{a1+a2+a5, -a5}, {-a5, a3+a4+a5}} *)
(* This is from the standard Symanzik form where Q = -[z^T Sigma z - 2J·z + C] *)

(* Let SA = a1+a2+a5, SB = a3+a4+a5 *)
(* Sigma = {{SA, -a5}, {-a5, SB}} *)

(* det(Sigma) = SA*SB - a5^2 = U *)
(* adj(Sigma) = {{SB, a5}, {a5, SA}} *)

SA = a1 + a2 + a5;
SB = a3 + a4 + a5;
Sigma = {{SA, -a5}, {-a5, SB}};
U = Det[Sigma];
adjSigma = {{SB, a5}, {a5, SA}};

Print["Sigma = ", Sigma // MatrixForm];
Print["U = det(Sigma) = ", U];
Print["adj(Sigma) = ", adjSigma // MatrixForm];
Print["Check: Sigma . adj(Sigma) = ", (Sigma . adjSigma) // MatrixForm, " = U * I ? ",
  Simplify[Sigma . adjSigma == U * IdentityMatrix[2]]];

(* ============================================================ *)
(* Step 4: Verify the F polynomial (X_i^2 cancellation) *)
(* ============================================================ *)

Print[""];
Print["Step 4: Verify F = C*U - J^T adj(Sigma) J, check X_i^2 cancellation"];

(* J = {a1*X1 + a2*X2, a3*X3 + a4*X4} *)
(* C = a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2 *)

(* We need to compute F = C*U - J^T adj(Sigma) J *)
(* and verify that all X_i^2 terms cancel *)

(* J^T adj(Sigma) J = (a1*X1+a2*X2)^T * SB * (a1*X1+a2*X2)
                     + 2*a5*(a1*X1+a2*X2)^T*(a3*X3+a4*X4)
                     + (a3*X3+a4*X4)^T * SA * (a3*X3+a4*X4) *)

(* Let's expand symbolically *)
(* Use symbolic X_i^2 and X_i·X_j *)
(* J^T adj(Sigma) J = *)
(*   SB * (a1^2*X1^2 + a2^2*X2^2 + 2*a1*a2*X1·X2) *)
(* + 2*a5 * (a1*a3*X1·X3 + a1*a4*X1·X4 + a2*a3*X2·X3 + a2*a4*X2·X4) *)
(* + SA * (a3^2*X3^2 + a4^2*X4^2 + 2*a3*a4*X3·X4) *)

(* C*U = (a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2) * (SA*SB - a5^2) *)
(* = (a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2) * SA*SB *)
(*   - (a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2) * a5^2 *)

(* Let's focus on X1^2 coefficient in F = C*U - J^T adj(Sigma) J *)
(* From C*U: a1 * SA*SB - a1 * a5^2 *)
(* From J^T adj(Sigma) J: a1^2 * SB *)
(* X1^2 coefficient: a1*SA*SB - a1*a5^2 - a1^2*SB *)
(* = a1*SB*(SA - a1) - a1*a5^2 *)
(* = a1*SB*(a2+a5) - a1*a5^2 *)
(* = a1*SB*a2 + a1*SB*a5 - a1*a5^2 *)
(* = a1*a2*SB + a1*a5*(SB - a5) *)
(* = a1*a2*SB + a1*a5*(a3+a4) *)

(* Hmm, this doesn't cancel! There are X1^2 terms remaining. *)
(* But the document says all X_i^2 terms cancel... *)

(* Wait, I think I made a sign error. Let me re-derive. *)
(* The standard form is: Q = sum(alpha_i * Delta_i) = -[z^T Sigma z - 2J·z + C] *)
(* Actually, Delta_i = -(X_i - z)^2, so Q = -sum(alpha_i * (X_i - z)^2) *)
(* This is Q = -[sum(alpha_i) * z^2 - 2*sum(alpha_i*X_i)·z + sum(alpha_i*X_i^2)] *)
(* ... for the single-z case. For 2 z's, it's more complex. *)

(* Let me re-derive from scratch. *)

(* Delta_1 = -(X1 - zA)^2 = -X1^2 + 2 X1·zA - zA^2 *)
(* Delta_2 = -(X2 - zA)^2 = -X2^2 + 2 X2·zA - zA^2 *)
(* Delta_5 = -(zA - zB)^2 = -zA^2 + 2 zA·zB - zB^2 *)
(* Delta_3 = -(X3 - zB)^2 = -X3^2 + 2 X3·zB - zB^2 *)
(* Delta_4 = -(X4 - zB)^2 = -X4^2 + 2 X4·zB - zB^2 *)

(* Q = a1*Delta_1 + a2*Delta_2 + a5*Delta_5 + a3*Delta_3 + a4*Delta_4 *)

(* = -(a1+a2+a5)*zA^2 - (a3+a4+a5)*zB^2 + 2*a5*zA·zB *)
(*   + 2*a1*X1·zA + 2*a2*X2·zA + 2*a3*X3·zB + 2*a4*X4·zB *)
(*   - a1*X1^2 - a2*X2^2 - a3*X3^2 - a4*X4^2 *)

(* This is: Q = -[z^T Sigma z - 2 J·z + C] *)
(* where: *)
(* z^T Sigma z = SA*zA^2 + SB*zB^2 - 2*a5*zA·zB *)
(* So Q = -[SA*zA^2 + SB*zB^2 - 2*a5*zA·zB] + 2*(a1*X1+a2*X2)·zA + 2*(a3*X3+a4*X4)·zB *)
(*         - (a1*X1^2+a2*X2^2+a3*X3^2+a4*X4^2) *)

(* Wait, the J term should be -2J·z in the bracket. *)
(* z^T Sigma z - 2J·z + C = SA*zA^2 + SB*zB^2 - 2*a5*zA·zB - 2*J1·zA - 2*J2·zB + C *)
(* = SA*zA^2 + SB*zB^2 - 2*a5*zA·zB - 2*(a1*X1+a2*X2)·zA - 2*(a3*X3+a4*X4)·zB + C *)

(* So -[z^T Sigma z - 2J·z + C] = -SA*zA^2 - SB*zB^2 + 2*a5*zA·zB + 2*(a1*X1+a2*X2)·zA + 2*(a3*X3+a4*X4)·zB - C *)

(* Comparing with Q: *)
(* Q = -(a1+a2+a5)*zA^2 - (a3+a4+a5)*zB^2 + 2*a5*zA·zB *)
(*   + 2*a1*X1·zA + 2*a2*X2·zA + 2*a3*X3·zB + 2*a4*X4·zB *)
(*   - a1*X1^2 - a2*X2^2 - a3*X3^2 - a4*X4^2 *)

(* So: SA = a1+a2+a5 ✓, SB = a3+a4+a5 ✓ *)
(* J1 = a1*X1 + a2*X2 ✓, J2 = a3*X3 + a4*X4 ✓ *)
(* C = a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2 ✓ *)

(* And the form is: Q = -[z^T Sigma z - 2J·z + C] *)

(* Now, the Gaussian integral: ∫ d^{2D}z Q^{-5a} *)
(* Complete the square: z^T Sigma z - 2J·z + C = (z - z0)^T Sigma (z - z0) + C - J^T Sigma^{-1} J *)
(* = (z - z0)^T Sigma (z - z0) + F/U *)
(* where F = C*U - J^T adj(Sigma) J *)

(* Now compute F: *)
(* C*U = (a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2) * (SA*SB - a5^2) *)

(* J^T adj(Sigma) J = SB*(a1*X1+a2*X2)^2 + 2*a5*(a1*X1+a2*X2)(a3*X3+a4*X4) + SA*(a3*X3+a4*X4)^2 *)

(* X1^2 coefficient in C*U: a1*SA*SB - a1*a5^2 *)
(* X1^2 coefficient in J^T adj(Sigma) J: a1^2*SB *)

(* X1^2 coefficient in F: a1*SA*SB - a1*a5^2 - a1^2*SB *)
(* = a1*SB*(SA - a1) - a1*a5^2 *)
(* SA - a1 = a2 + a5 *)
(* = a1*SB*(a2+a5) - a1*a5^2 *)
(* = a1*a2*SB + a1*a5*SB - a1*a5^2 *)
(* = a1*a2*SB + a1*a5*(SB - a5) *)
(* SB - a5 = a3 + a4 *)
(* = a1*a2*SB + a1*a5*(a3+a4) *)

(* This is NOT zero! X1^2 doesn't cancel in F. *)

(* But the document says F = -B_5, and B_5 has no X_i^2 terms. *)
(* This means either: *)
(* 1. My derivation is wrong *)
(* 2. The document is wrong *)
(* 3. There's an extra step I'm missing *)

(* Actually, wait. Let me re-read the document. It says: *)
(* "展开后所有 X_i^2 项消去" *)
(* This is a claim about F, not about C*U - J^T adj(Sigma) J. *)

(* Let me re-derive F more carefully. *)
(* After completing the square: *)
(* Q = -[(z - z0)^T Sigma (z - z0) + F/U] *)
(* = -(z - z0)^T Sigma (z - z0) - F/U *)

(* So Q = - (positive quadratic form) - F/U *)

(* The Gaussian integral: ∫ d^{2D}z Q^{-5a} = ∫ d^{2D}z [-(z-z0)^T Sigma (z-z0) - F/U]^{-5a} *)
(* = (-1)^{-5a} ∫ d^{2D}z [(z-z0)^T Sigma (z-z0) + F/U]^{-5a} *)

(* Hmm, this is getting complicated. Let me think about this differently. *)

(* Actually, looking at the symanzik-parametrization.md document: *)
(* "Q = -z^T(Sigma ⊗ I_D)z + 2J·z - C" *)
(* Wait, that's different from what I wrote. Let me re-check. *)

(* The symanzik doc says: Q = z^T(Sigma ⊗ I_D)z - 2J·z + C *)
(* But Delta_i = -(X_i - z)^2, so Q = sum(alpha_i * (-(X_i - z)^2)) *)
(* = -sum(alpha_i * (X_i^2 - 2X_i·z + z^2)) *)
(* = -sum(alpha_i)*z^2 + 2*sum(alpha_i*X_i)·z - sum(alpha_i*X_i^2) *)
(* = z^T*(-sum(alpha_i)*I)z + 2*J·z - C *)
(* So Q = z^T(-I)z + 2J·z - C for the 1-z case. *)
(* This is NOT z^T Sigma z - 2J·z + C. The signs are different. *)

(* Let me check the symanzik doc again more carefully. *)
(* "Q = z^T(Sigma ⊗ I_D)z - 2J·z + C" *)
(* With Sigma_11 = a1+a2+a5 for the 5-point case. *)
(* And Q = a1*Delta_1 + ... + a5*Delta_5 *)

(* Delta_1 = -(X1 - zA)^2 = -(X1^2 - 2X1·zA + zA^2) = -X1^2 + 2X1·zA - zA^2 *)

(* Q = a1*(-X1^2 + 2X1·zA - zA^2) + a2*(-X2^2 + 2X2·zA - zA^2) + a5*(-zA^2 + 2zA·zB - zB^2) *)
(*     + a3*(-X3^2 + 2X3·zB - zB^2) + a4*(-X4^2 + 2X4·zB - zB^2) *)

(* = -(a1+a2+a5)*zA^2 - (a3+a4+a5)*zB^2 + 2*a5*zA·zB *)
(*   + 2a1*X1·zA + 2a2*X2·zA + 2a3*X3·zB + 2a4*X4·zB *)
(*   - (a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2) *)

(* For the form Q = z^T Sigma z - 2J·z + C: *)
(* z^T Sigma z = Sigma_11*zA^2 + 2*Sigma_12*zA·zB + Sigma_22*zB^2 *)
(* -2J·z = -2*J1·zA - 2*J2·zB *)

(* Matching: Sigma_11 = -(a1+a2+a5) = -SA *)
(* 2*Sigma_12 = 2*a5 → Sigma_12 = a5 *)
(* Sigma_22 = -(a3+a4+a5) = -SB *)
(* J1 = -(a1*X1 + a2*X2) *)
(* J2 = -(a3*X3 + a4*X4) *)
(* C = -(a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2) *)

(* So Sigma = {{-SA, a5}, {a5, -SB}} *)

(* But the document says Sigma = {{SA, -a5}, {-a5, SB}}. *)
(* That's Sigma = -{{-SA, a5}, {a5, -SB}}. The document's Sigma is the NEGATIVE of mine! *)

(* Actually, I think the convention in the symanzik doc is: *)
(* Q = -[z^T Sigma z - 2J·z + C]  (note the outer minus sign) *)
(* So Q = -z^T Sigma z + 2J·z - C *)
(* = z^T (-Sigma) z + 2J·z - C *)

(* With Sigma (document) = {{SA, -a5}, {-a5, SB}}: *)
(* -Sigma = {{-SA, a5}, {a5, -SB}} *)
(* This matches my coefficients! ✓ *)

(* So the convention is Q = -[z^T Sigma z - 2J·z + C] with Sigma = {{SA, -a5}, {-a5, SB}}. *)

(* After completing the square: *)
(* z^T Sigma z - 2J·z + C = (z - z0)^T Sigma (z - z0) + C - J^T Sigma^{-1} J *)
(* = (z - z0)^T Sigma (z - z0) + F/U *)

(* Q = -[(z - z0)^T Sigma (z - z0) + F/U] *)

(* Gaussian integral: ∫ d^{2D}z Q^{-5a} = ∫ d^{2D}z [-(z-z0)^T Sigma (z-z0) - F/U]^{-5a} *)
(* = (-1)^{-5a} ∫ d^{2D}z [(z-z0)^T Sigma (z-z0) + F/U]^{-5a} *)

(* Actually, (-1)^{-5a} = (-1)^{-5(1-ε)} = (-1)^{-5+5ε} = (-1)^{-5} * (-1)^{5ε} *)
(* (-1)^{-5} = -1 *)
(* This is where the minus sign in the prefactor comes from! *)

(* Now, the integral ∫ d^{2D}z [(z-z0)^T Sigma (z-z0) + F/U]^{-5a} *)
(* Shift z → z+z0, then the integral is standard 2D Gaussian. *)
(* = (iπ^D)^2? No, wait. *)

(* Let me use the standard formula: *)
(* ∫ d^{2D}z [z^T (Sigma ⊗ I_D) z + F/U]^{-5a} *)
(* = π^D / Gamma(5a) * Gamma(5a-D) * U^{-D/2} * (F/U)^{-(5a-D)} *)

(* Hmm, I keep getting confused. Let me just verify the key claim: *)
(* does F = C*U - J^T adj(Sigma) J have X_i^2 terms? *)

(* I already computed the X1^2 coefficient in F: *)
(* a1*a2*SB + a1*a5*(a3+a4) *)
(* This is NOT zero. *)

(* But wait — F is NOT B_5. The document says B_5 = -F. *)
(* So B_5 would have X_i^2 terms too, which contradicts the document's B_5 formula. *)

(* Let me double-check by computing F explicitly with a symbolic algebra approach. *)

(* Actually, let me just verify the 3-point case first to make sure my method is correct. *)
Print[""];
Print["Step 4a: First verify 3-point case (known result)"];

(* 3-point: 1 z, 3 Delta_i = -(X_i - z)^2 *)
(* Sigma = a1 + a2 + a3 = 1 (with delta constraint) *)
(* J = a1*X1 + a2*X2 + a3*X3 *)
(* C = a1*X1^2 + a2*X2^2 + a3*X3^2 *)
(* U = 1 *)
(* F = C*1 - J^2 *)

(* C = a1*X1^2 + a2*X2^2 + a3*X3^2 *)
(* J^2 = a1^2*X1^2 + a2^2*X2^2 + a3^2*X3^2 + 2a1*a2*X1·X2 + 2a2*a3*X2·X3 + 2a3*a1*X3·X1 *)

(* F = (a1 - a1^2)*X1^2 + (a2 - a2^2)*X2^2 + (a3 - a3^2)*X3^2 *)
(*     - 2a1*a2*X1·X2 - 2a2*a3*X2·X3 - 2a3*a1*X3·X1 *)

(* Using a1 + a2 + a3 = 1: *)
(* a1 - a1^2 = a1*(1-a1) = a1*(a2+a3) *)
(* So X1^2 term: a1*a2*X1^2 + a1*a3*X1^2 *)

(* F = a1*a2*(X1^2 + X2^2 - 2X1·X2) + a2*a3*(X2^2 + X3^2 - 2X2·X3) + a3*a1*(X3^2 + X1^2 - 2X3·X1) *)
(*   = a1*a2*(X1-X2)^2 + a2*a3*(X2-X3)^2 + a3*a1*(X3-X1)^2 *)
(*   = -a1*a2*Delta_12 - a2*a3*Delta_23 - a3*a1*Delta_31 *)

(* B_3 = -F = a1*a2*Delta_12 + a2*a3*Delta_23 + a3*a1*Delta_31 ✓ *)
(* All X_i^2 terms DID cancel! *)

Print["3-point: F = C - J^2 = -B_3. All X_i^2 cancel ✓ (verified algebraically)"];

(* Now do the 5-point case properly *)
Print[""];
Print["Step 4b: 5-point F computation"];

(* U = SA*SB - a5^2 *)
(* C = a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2 *)
(* J = {a1*X1 + a2*X2, a3*X3 + a4*X4} *)
(* adj(Sigma) = {{SB, a5}, {a5, SA}} *)

(* J^T adj(Sigma) J = SB*(a1*X1+a2*X2)^2 + 2*a5*(a1*X1+a2*X2)(a3*X3+a4*X4) + SA*(a3*X3+a4*X4)^2 *)

(* X1^2 term in F: *)
(* From C*U: a1*SA*SB - a1*a5^2 *)
(* From J^T adj(Sigma) J: a1^2*SB *)
(* X1^2 coefficient: a1*SA*SB - a1*a5^2 - a1^2*SB *)
(* = a1*SB*(SA - a1) - a1*a5^2 *)
(* = a1*SB*(a2+a5) - a1*a5^2 *)
(* = a1*a2*SB + a1*a5*SB - a1*a5^2 *)
(* = a1*a2*SB + a1*a5*(a3+a4+a5) - a1*a5^2 *)
(* = a1*a2*SB + a1*a5*(a3+a4) *)

(* This is NOT zero! X1^2 terms remain in F. *)
(* But the document claims B_5 = -F has no X_i^2 terms. *)

(* Wait, am I computing F with the right sign convention? *)
(* F = C*U - J^T adj(Sigma) J. *)
(* After completing the square: Q = -[(z-z0)^T Sigma (z-z0) + F/U] *)
(* So Q = -(z-z0)^T Sigma (z-z0) - F/U *)

(* The Gaussian integral: *)
(* ∫ d^{2D}z Q^{-5a} = ∫ d^{2D}z [-(z-z0)^T Sigma (z-z0) - F/U]^{-5a} *)
(* Let z' = z - z0, then: *)
(* = ∫ d^{2D}z' [-(z'^T Sigma z') - F/U]^{-5a} *)
(* = (-1)^{-5a} ∫ d^{2D}z' [(z'^T Sigma z') + F/U]^{-5a} *)

(* For the integral to converge, we need (z'^T Sigma z') + F/U > 0. *)
(* Sigma is not positive definite in general... *)

(* Actually, I think the issue is that the standard Symanzik formula assumes *)
(* Q = z^T Sigma z - 2J·z + C (without the outer minus sign), and the Delta_i *)
(* are defined as +(X_i - z)^2, not -(X_i - z)^2. *)

(* Let me check the standard convention. In the Symanzik parametrization, the *)
(* propagators are 1/(-k^2 + i0). After Fourier transform, we get 1/(-(x)^2 + i0)^a. *)
(* So Delta_i = -(X_i - z)^2 + i0. Then Q = sum(alpha_i * Delta_i). *)

(* In the Gaussian integral, we have Q^{-5a} = [sum(alpha_i * Delta_i)]^{-5a}. *)
(* The standard formula is for ∫ d^Dz [z^T A z + 2b·z + c]^{-N} *)

(* Let me use a different approach: directly verify the 5-point identity *)
(* by computing the integral for a simple case where all X_i = 0. *)
(* Actually, that's degenerate. Let me use a numeric verification approach. *)

(* Actually, let me try a different approach. Let me compute F directly using *)
(* symbolic algebra, treating X_i^2 and X_i·X_j as independent symbols. *)

Print[""];
Print["Step 4c: Direct symbolic computation of F"];

(* Let's use Mathematica's symbolic capabilities *)
(* Treat X1sq = X1^2, X2sq = X2^2, etc. as symbols *)
(* And X1X2 = X1·X2, etc. *)

(* C = a1*X1sq + a2*X2sq + a3*X3sq + a4*X4sq *)
(* J1 = a1*X1 + a2*X2, but we need J1^2 = (a1*X1 + a2*X2)^2 = a1^2*X1sq + a2^2*X2sq + 2*a1*a2*X1X2 *)
(* Similarly for J2^2 and J1·J2 *)

(* Let me compute F = C*U - J^T adj(Sigma) J symbolically *)

(* C*U = (a1*X1sq + a2*X2sq + a3*X3sq + a4*X4sq) * (SA*SB - a5^2) *)

(* J^T adj(Sigma) J = *)
(* SB * (a1^2*X1sq + a2^2*X2sq + 2*a1*a2*X1X2) *)
(* + 2*a5 * (a1*a3*X1X3 + a1*a4*X1X4 + a2*a3*X2X3 + a2*a4*X2X4) *)
(* + SA * (a3^2*X3sq + a4^2*X4sq + 2*a3*a4*X3X4) *)

(* X1sq coefficient in F: *)
(* From C*U: a1*(SA*SB - a5^2) *)
(* From J^T adj(Sigma) J: a1^2*SB *)
(* Total: a1*SA*SB - a1*a5^2 - a1^2*SB = a1*SB*(SA - a1) - a1*a5^2 *)
(* = a1*SB*(a2+a5) - a1*a5^2 = a1*a2*SB + a1*a5*SB - a1*a5^2 *)
(* = a1*a2*SB + a1*a5*(SB - a5) = a1*a2*SB + a1*a5*(a3+a4) *)

(* This is NOT zero. X1^2 does NOT cancel from F. *)

(* Hmm, but the symanzik-parametrization.md doc clearly says: *)
(* "Python 展开（scripts/compute_F_polynomial.py 验证：所有 X_i^2 项消去，F 纯由 Δ_ij 表达）" *)

(* Let me check if that Python script exists and what it computes. *)
(* But first, let me re-check my algebra. *)

(* Actually, wait. Maybe I'm confusing F and B_5. Let me re-read: *)
(* B_5 = -F, and the document says B_5 has no X_i^2 terms. *)
(* But my computation shows F (and hence B_5) has X_i^2 terms. *)

(* UNLESS... the definition of F in the document is different from what I'm computing. *)
(* Let me re-read the symanzik doc: *)
(* "F = C·U - J^T adj(Sigma) J" *)
(* "B_5 = -F" *)

(* But my computation shows X1^2 coefficient = a1*a2*SB + a1*a5*(a3+a4) ≠ 0 *)

(* Actually wait, let me reconsider. Perhaps the F polynomial is defined *)
(* differently. In the standard Symanzik formalism: *)
(* Q = z^T (Sigma ⊗ I_D) z - 2 J·z + C *)
(* Complete square: Q = (z - z0)^T (Sigma ⊗ I_D) (z - z0) - F/U *)
(* where F = J^T adj(Sigma) J - C*U  (note the sign difference!) *)

(* Let me check this. If Q = z^T Sigma z - 2J·z + C: *)
(* z0 = Sigma^{-1} J *)
(* Q = (z - z0)^T Sigma (z - z0) + C - J^T Sigma^{-1} J *)
(* = (z - z0)^T Sigma (z - z0) + C - J^T (adj(Sigma)/U) J *)
(* = (z - z0)^T Sigma (z - z0) + (C*U - J^T adj(Sigma) J)/U *)

(* So F = C*U - J^T adj(Sigma) J. *)
(* And Q = (z - z0)^T Sigma (z - z0) + F/U. *)

(* Then the Gaussian integral is: *)
(* ∫ d^{2D}z Q^{-5a} = ∫ d^{2D}z [(z - z0)^T Sigma (z - z0) + F/U]^{-5a} *)

(* For this to work, we need Q = (positive) + F/U. But our Q = sum(alpha_i * Delta_i) *)
(* and Delta_i = -(X_i - z)^2, so Q is negative definite in z. *)
(* This means Sigma should be negative definite, and F/U should be negative too. *)

(* Actually, in the standard convention, Q = z^T Sigma z - 2J·z + C *)
(* and the propagators are taken as 1/(-k^2), so after FT we get positive Delta. *)
(* But our convention uses Delta = -(X_i - z)^2, which gives a negative Q. *)

(* The resolution: the standard Symanzik formula is for Q = z^T A z + 2b·z + c *)
(* (note the PLUS sign for the linear term). Our Q has: *)
(* Q = -(SA*zA^2 + SB*zB^2 - 2*a5*zA·zB) + 2*(a1*X1+a2*X2)·zA + 2*(a3*X3+a4*X4)·zB *)
(*     - (a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2) *)

(* = -z^T Sigma z + 2J·z - C *)
(* where Sigma = {{SA, -a5}, {-a5, SB}}, J = (a1*X1+a2*X2, a3*X3+a4*X4), C = sum(a_i*X_i^2) *)

(* Complete square: -z^T Sigma z + 2J·z - C = -(z^T Sigma z - 2J·z + C) *)
(* = -[(z - z0)^T Sigma (z - z0) + (C*U - J^T adj(Sigma) J)/U] *)
(* = -(z - z0)^T Sigma (z - z0) - F/U *)

(* So Q = -(z - z0)^T Sigma (z - z0) - F/U *)
(* and Q^{-5a} = (-1)^{-5a} [(z - z0)^T Sigma (z - z0) + F/U]^{-5a} *)

(* Now the Gaussian integral formula applies to [(z - z0)^T Sigma (z - z0) + F/U]^{-5a}. *)
(* The result is proportional to (F/U)^{-(5a-D)} * U^{-D/2} = (F)^{-(5a-D)} * U^{5a-3D/2}. *)

(* Wait, I need to check the exponents more carefully. *)
(* The standard formula: ∫ d^{LD}z [z^T A z + m^2]^{-N} = π^{LD/2} Γ(N-LD/2)/Γ(N) * det(A)^{-D/2} * (m^2)^{-(N-LD/2)} *)

(* Here A = Sigma ⊗ I_D, det(A) = U^D, N = 5a, LD = 2D, m^2 = F/U. *)
(* Result: π^D Γ(5a-D)/Γ(5a) * U^{-D/2} * (F/U)^{-(5a-D)} *)
(* = π^D Γ(5a-D)/Γ(5a) * U^{-D/2 - (5a-D)} * F^{-(5a-D)}? No... *)

(* (F/U)^{-(5a-D)} = F^{-(5a-D)} * U^{5a-D} *)
(* Total U exponent: -D/2 + 5a - D = 5a - 3D/2 = 5(1-ε) - 3(4-2ε)/2 = 5-5ε-6+3ε = -1-2ε *)
(* Wait, that's U^{-(1+2ε)}. And the document has U^{-(1+2ε)}. ✓ *)

(* F exponent: -(5a-D) = -(5-5ε-4+2ε) = -(1-3ε) *)
(* Wait, 5a - D = 5(1-ε) - (4-2ε) = 1 - 3ε *)
(* So F^{-(1-3ε)}. And the document has B_5^{-(1-3ε)}. *)
(* Since B_5 = -F, B_5^{-(1-3ε)} = (-F)^{-(1-3ε)} = (-1)^{-(1-3ε)} F^{-(1-3ε)}. *)
(* Hmm, there's an extra (-1) factor that needs to be tracked. *)

(* OK, this is getting very involved. Let me just verify the key algebraic claim: *)
(* does F = C*U - J^T adj(Sigma) J have X_i^2 terms? *)

(* My computation says YES. But the document says NO. *)
(* Let me write a concrete numeric check to resolve this. *)

Print[""];
Print["Step 4d: Numeric check with random values"];

(* Random values for alpha_i (satisfying sum = 1) *)
r[] := RandomReal[{0,1}];
Module[{a1v, a2v, a3v, a4v, a5v, sum},
  a1v = r[]; a2v = r[]; a3v = r[]; a4v = r[]; a5v = r[];
  sum = a1v + a2v + a3v + a4v + a5v;
  a1v /= sum; a2v /= sum; a3v /= sum; a4v /= sum; a5v /= sum;

  SAv = a1v + a2v + a5v;
  SBv = a3v + a4v + a5v;
  Uv = SAv*SBv - a5v^2;

  (* Random X_i^2 and X_i·X_j values *)
  x1sq = r[]; x2sq = r[]; x3sq = r[]; x4sq = r[];
  x1x2 = r[]; x1x3 = r[]; x1x4 = r[]; x2x3 = r[]; x2x4 = r[]; x3x4 = r[];

  (* C * U *)
  CU = (a1v*x1sq + a2v*x2sq + a3v*x3sq + a4v*x4sq) * Uv;

  (* J^T adj(Sigma) J *)
  JTadjJ = SBv*(a1v^2*x1sq + a2v^2*x2sq + 2*a1v*a2v*x1x2) +
           2*a5v*(a1v*a3v*x1x3 + a1v*a4v*x1x4 + a2v*a3v*x2x3 + a2v*a4v*x2x4) +
           SAv*(a3v^2*x3sq + a4v^2*x4sq + 2*a3v*a4v*x3x4);

  Fv = CU - JTadjJ;

  (* Compute expected B_5 from the document's formula *)
  (* Delta_ij = -(X_i - X_j)^2 + i0 = -X_i^2 - X_j^2 + 2*X_i·X_j + i0 *)
  (* But for the X_i^2 check, we just need to verify the structure *)
  (* B_5 = a1*a2*SB*Delta_12 + a3*a4*SA*Delta_34 + a5*(a1*a3*Delta_13 + ...) *)

  (* Let's compute B_5 in terms of X_i^2 and X_i·X_j *)
  delta12 = -(x1sq + x2sq - 2*x1x2);  (* Delta_12 = -(X1-X2)^2 *)
  delta34 = -(x3sq + x4sq - 2*x3x4);
  delta13 = -(x1sq + x3sq - 2*x1x3);
  delta14 = -(x1sq + x4sq - 2*x1x4);
  delta23 = -(x2sq + x3sq - 2*x2x3);
  delta24 = -(x2sq + x4sq - 2*x2x4);

  B5v = a1v*a2v*SBv*delta12 + a3v*a4v*SAv*delta34 +
        a5v*(a1v*a3v*delta13 + a1v*a4v*delta14 + a2v*a3v*delta23 + a2v*a4v*delta24);

  Print["F = ", Fv];
  Print["-B_5 = ", -B5v];
  Print["F + B_5 = ", Fv + B5v];
  Print["F = -B_5 ? ", If[Abs[Fv + B5v] < 10^-10, "✓ YES", "✗ NO"]];
];

(* ============================================================ *)
Print[""];
Print["===== VERIFICATION SUMMARY ====="];
Print["Sigma structure: verified algebraically"];
Print["U = det(Sigma): verified algebraically"];
Print["adj(Sigma): verified algebraically"];
Print["F = C*U - J^T adj(Sigma) J: pending numeric check above"];
Print["Gaussian integral formula: standard result, prefactor to be verified separately"];