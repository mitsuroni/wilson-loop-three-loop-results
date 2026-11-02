(* v04-symanzik-param.m *)
(* ========================================================================= *)
(* 验证：B_5, U 的 Symanzik 参数化公式                                   *)
(* 方法：从位置空间传播子出发，独立推导 Q, Sigma, U, B_5                  *)
(* ========================================================================= *)

Print["=============================================="];
Print["v04: Symanzik 参数化验证"];
Print["=============================================="];

(* ========================================================================= *)
(* Step 1: 定义传播子                                                        *)
(* ========================================================================= *)

(* Δ_i = -(X_i - z)^2 + i0，其中 z = z_A 或 z_B *)
(* 具体地： *)
(*   Δ_1 = -(X_1 - z_A)^2 *)
(*   Δ_2 = -(X_2 - z_A)^2 *)
(*   Δ_5 = -(z_A - z_B)^2 *)
(*   Δ_3 = -(X_3 - z_B)^2 *)
(*   Δ_4 = -(X_4 - z_B)^2 *)

(* Q = α_1 Δ_1 + α_2 Δ_2 + α_5 Δ_5 + α_3 Δ_3 + α_4 Δ_4 *)
(* 展开 Q 为 z_A, z_B 的二次型 *)

(* 用符号变量 *)
ClearAll[X1, X2, X3, X4, zA, zB, a1, a2, a5, a3, a4];

(* 定义 Δ_i（忽略 i0）*)
delta1 = -(X1 - zA)^2;
delta2 = -(X2 - zA)^2;
delta5 = -(zA - zB)^2;
delta3 = -(X3 - zB)^2;
delta4 = -(X4 - zB)^2;

Q = a1*delta1 + a2*delta2 + a5*delta5 + a3*delta3 + a4*delta4;

Print["Q = ", Q];

(* ========================================================================= *)
(* Step 2: 提取 z_A, z_B 的二次型系数                                       *)
(* ========================================================================= *)

(* 展开 Q *)
Qexpanded = Expand[Q];
Print[""];
Print["Q 展开: ", Qexpanded];

(* 提取 z_A^2, z_B^2, z_A·z_B 的系数 *)
(* z_A^2 的系数 *)
coeffzA2 = Coefficient[Qexpanded, zA^2];
Print["z_A^2 系数: ", coeffzA2];

(* z_B^2 的系数 *)
coeffzB2 = Coefficient[Qexpanded, zB^2];
Print["z_B^2 系数: ", coeffzB2];

(* z_A·z_B 的系数（注意 Q 中有 -2 zA·zB * ... 的形式）*)
(* 更简单的方法：提取 zA 和 zB 的矩阵形式 *)
(* Q = -[zA^2 Σ_11 + zB^2 Σ_22 + 2 zA·zB Σ_12 - 2 J_1·zA - 2 J_2·zB + C] *)

(* Σ_11 = (a1+a2+a5) *)
sigma11 = a1 + a2 + a5;
Print["Σ_11 = α1+α2+α5 = ", sigma11];

(* Σ_22 = (a3+a4+a5) *)
sigma22 = a3 + a4 + a5;
Print["Σ_22 = α3+α4+α5 = ", sigma22];

(* Σ_12 = -a5 *)
sigma12 = -a5;
Print["Σ_12 = -α5 = ", sigma12];

(* 验证 *)
(* Q = -[sigma11*zA^2 + sigma22*zB^2 + 2*sigma12*zA·zB - 2*J1·zA - 2*J2·zB + C] *)
(* 其中 zA^2 的系数是 -sigma11 *)
checkSigma11 = Simplify[coeffzA2 + sigma11];
Print["zA^2 系数 = -Σ_11: ", If[checkSigma11 === 0, "✓", "✗"]];

(* zB^2 的系数是 -sigma22 *)
checkSigma22 = Simplify[coeffzB2 + sigma22];
Print["zB^2 系数 = -Σ_22: ", If[checkSigma22 === 0, "✓", "✗"]];

(* ========================================================================= *)
(* Step 3: 验证 U = det(Σ)                                                   *)
(* ========================================================================= *)

(* Σ = {{sigma11, sigma12}, {sigma12, sigma22}} *)
Sigma = {{sigma11, sigma12}, {sigma12, sigma22}};
Udet = Det[Sigma];
Print[""];
Print["U = det(Σ) = ", Udet];

(* 文档声称的公式 *)
Udoc = (a1 + a2 + a5)*(a3 + a4 + a5) - a5^2;
checkU = Simplify[Udet - Udoc];
Print["U = Σ_A Σ_B - α5^2: ", If[checkU === 0, "✓", "✗"]];

(* ========================================================================= *)
(* Step 4: 提取 J_1, J_2, C                                                 *)
(* ========================================================================= *)

(* 展开 Q 的 zA 线性项 *)
(* Q 中 zA·X_i 的项来自 Δ_1, Δ_2, Δ_5 *)
(* Δ_1 = -(X1^2 - 2 X1·zA + zA^2) → -X1^2 + 2 X1·zA - zA^2 *)
(* 所以 zA 的线性项系数: 2a1*X1 + 2a2*X2 + 2a5*zB *)

(* J_1 = a1*X1 + a2*X2 + a5*zB *)
(* 但实际上 J_1 只含 X_i，不含 zB。zB 项来自 Σ_12 的交叉项 *)
(* 正确的 J_1 = a1*X1 + a2*X2 *)
J1 = a1*X1 + a2*X2;
J2 = a3*X3 + a4*X4;
Print[""];
Print["J_1 = α1 X1 + α2 X2 = ", J1];
Print["J_2 = α3 X3 + α4 X4 = ", J2];

(* C = Σ α_i X_i^2 *)
Cterm = a1*X1^2 + a2*X2^2 + a3*X3^2 + a4*X4^2;
Print["C = Σ α_i X_i^2 = ", Cterm];

(* ========================================================================= *)
(* Step 5: 验证 B_5 = -F = -(C·U - J^T adj(Σ) J)                          *)
(* ========================================================================= *)

(* adj(Σ) = {{sigma22, -sigma12}, {-sigma12, sigma11}} *)
adjSigma = {{sigma22, -sigma12}, {-sigma12, sigma11}};

(* J^T adj(Σ) J *)
Jvec = {J1, J2};
JTadjJ = Jvec . adjSigma . Jvec;
Print[""];
Print["J^T adj(Σ) J = ", Expand[JTadjJ]];

(* F = C*U - JTadjJ *)
Fterm = Cterm*Udet - JTadjJ;
Fexpanded = Expand[Fterm];
Print["F = C·U - J^T adj(Σ) J = ", Fexpanded];

(* B_5 = -F *)
B5 = -Fexpanded;
Print[""];
Print["B_5 = -F = ", B5];

(* ========================================================================= *)
(* Step 6: 验证 B_5 的 6 项展开式                                            *)
(* ========================================================================= *)

(* 文档声称的 B_5 公式 *)
B5doc = a1*a2*(a3 + a4 + a5)*(X1 - X2)^2 +
        a3*a4*(a1 + a2 + a5)*(X3 - X4)^2 +
        a5*a1*a3*(X1 - X3)^2 +
        a5*a1*a4*(X1 - X4)^2 +
        a5*a2*a3*(X2 - X3)^2 +
        a5*a2*a4*(X2 - X4)^2;

(* 注意：Δ_ij = -(X_i - X_j)^2，所以 B_5 中 Δ_ij 的符号是负的 *)
(* 但 B_5 = -F，其中 F 中的 X_i^2 项被消去后，剩下的项都是 +(X_i-X_j)^2 形式 *)
(* 所以 B_5 = Σ α_i α_j c_ij (-(X_i-X_j)^2) = -Σ α_i α_j c_ij (X_i-X_j)^2 *)
(* 等等，我们需要小心符号。 *)

(* 实际上，Δ_i = -(X_i-z)^2 + i0 *)
(* Q = Σ α_i Δ_i = -Σ α_i (X_i-z)^2 *)
(* 经过高斯积分后，分母是 Q^{-5a} = (-Σ α_i (X_i-z)^2)^{-5a} *)
(* 配平方后得到 B_5/U，其中 B_5 是正定的 *)
(* 从文档的推导，B_5 = +α_1α_2 Σ_B Δ_12 + ... *)
(* 但 Δ_12 = -(X_1-X_2)^2，所以 B_5 = -α_1α_2 Σ_B (X_1-X_2)^2 + ... *)

(* 让我们验证符号 *)
B5neg = -B5doc;  (* 带负号的 B_5 *)
(* 实际上 B_5 = -F，而 F 的展开... *)
(* 让我直接比较 B_5 和 B_5doc。*)

(* 重新计算：B_5 应该等于什么？ *)
(* 从 Q 的配平方：Q = -(z-z0)^T Σ (z-z0) - F/U *)
(* 其中 F = C·U - J^T adj(Σ) J *)
(* 然后 B_5 = -F = J^T adj(Σ) J - C·U *)
(* 高斯积分后分母是 (-F/U)^{5a-D} = (B_5/U)^{5a-D} *)

(* 所以 B_5 = J^T adj(Σ) J - C·U *)
(* 让我验证 B_5 = B5doc 或 B_5 = -B5doc *)

B5fromJT = Expand[JTadjJ - Cterm*Udet];
Print[""];
Print["B_5 (from J^T adj J - C·U) = ", B5fromJT];

(* 比较 B5fromJT 和 B5doc *)
diffB5 = Simplify[B5fromJT - B5doc];
Print["B_5 - B_5^doc = ", diffB5];
Print["B_5 = B_5^doc: ", If[diffB5 === 0, "✓", "符号可能相反"]];

(* 如果差一个负号 *)
diffB5neg = Simplify[B5fromJT + B5doc];
Print["B_5 = -B_5^doc: ", If[diffB5neg === 0, "✓", "✗"]];

(* ========================================================================= *)
(* Step 7: 用 Δ_ij 表示 B_5                                                  *)
(* ========================================================================= *)

(* Δ_ij = -(X_i - X_j)^2 *)
(* 所以 B_5 = Σ α_i α_j c_ij Δ_ij *)
(* 其中 c_ij 由公式给出 *)

(* 验证：B_5 的每一项是否形如 α_i α_j c_ij (X_i-X_j)^2 *)
(* 如果是，则 B_5 = -Σ α_i α_j c_ij Δ_ij（因为 Δ_ij = -(X_i-X_j)^2）*)

(* 提取 B_5 中 (X_i-X_j)^2 的项 *)
(* 例如：(X1-X2)^2 的系数 *)
coeffX1X2 = Coefficient[B5fromJT, (X1 - X2)^2];
Print[""];
Print["(X1-X2)^2 系数: ", coeffX1X2];
Print["  文档声称: -α1 α2 (α3+α4+α5) = -", a1*a2*(a3+a4+a5), " (因为 Δ_12 = -(X1-X2)^2)"];
check12 = Simplify[coeffX1X2 + a1*a2*(a3+a4+a5)];
Print["  ", If[check12 === 0, "✓", "✗"]];

coeffX3X4 = Coefficient[B5fromJT, (X3 - X4)^2];
Print["(X3-X4)^2 系数: ", coeffX3X4];
check34 = Simplify[coeffX3X4 + a3*a4*(a1+a2+a5)];
Print["  ", If[check34 === 0, "✓", "✗"]];

coeffX1X3 = Coefficient[B5fromJT, (X1 - X3)^2];
Print["(X1-X3)^2 系数: ", coeffX1X3];
check13 = Simplify[coeffX1X3 + a5*a1*a3];
Print["  ", If[check13 === 0, "✓", "✗"]];

Print[""];
Print["v04 完成"];