(* v07-orbit25-explicit.m *)
(* ========================================================================= *)
(* orbit 25: (1,2,3,4) 非退化 — 完整 C × J 显式计算                     *)
(* ========================================================================= *)

Print["=============================================="];
Print["v07: orbit 25 (1,2,3,4) 完整 C x J"];
Print["=============================================="];

(* ========================================================================= *)
(* Step 1: 运动学与几何                                                      *)
(* ========================================================================= *)

(* N=4 运动学 *)
ClearAll[z, vRef];
(* 用列表表示 v_i·v_j / vRef *)
vdotTable = <|
  {1,2} -> z, {3,4} -> z,
  {2,3} -> 1, {1,4} -> 1,
  {1,3} -> -(1+z), {2,4} -> -(1+z)
|>;
vdot[i_, j_] := If[i==j, 0, vdotTable[Sort[{i,j}]]] * vRef;

(* 边分布: a=1, b=2, c=3, d=4 *)
a=1; b=2; c=3; d=4;

(* 位置: X_i = V_i + s_i v_i *)
(* V_1 = 0, V_2 = v_1, V_3 = v_1+v_2, V_4 = v_1+v_2+v_3 *)
(* 用符号 s1, s2, s3, s4 代替 s_{1,1}, s_{2,1}, s_{3,1}, s_{4,1} *)

(* 计算 v_p·(X_i - X_j) / vRef *)
(* X_i - X_j = (V_i - V_j) + s_i v_i - s_j v_j *)
(* v_p·(X_i - X_j) = v_p·(V_i - V_j) + s_i v_p·v_i - s_j v_p·v_j *)

ClearAll[s1, s2, s3, s4];
vDotXDifference[p_, i_, j_] := Module[
  {vi, vj, si, sj, Vdiff},
  vi = Switch[i, 1,a, 2,b, 3,c, 4,d];
  vj = Switch[j, 1,a, 2,b, 3,c, 4,d];
  si = Switch[i, 1,s1, 2,s2, 3,s3, 4,s4];
  sj = Switch[j, 1,s1, 2,s2, 3,s3, 4,s4];
  (* V_i - V_j *)
  Vdiff = vdot[p, vi]*si/vRef - vdot[p, vj]*sj/vRef;
  (* 加上 v_p·(V_i-V_j)/vRef *)
  (* V_1=0, V_2=v_1, V_3=v_1+v_2, V_4=v_1+v_2+v_3 *)
  Module[{Vi, Vj},
    Vi = Switch[i, 1, {0}, 2, {1}, 3, {1,2}, 4, {1,2,3}];
    Vj = Switch[j, 1, {0}, 2, {1}, 3, {1,2}, 4, {1,2,3}];
    (* V_i - V_j = 从 V_j 到 V_i 的边向量之和 *)
    (* 简单方法：V_k = sum_{m=1}^{k-1} v_m *)
    (* V_i - V_j = sum_{m=j}^{i-1} v_m (if i>j) or -sum_{m=i}^{j-1} v_m (if i<j) *)
    Module[{lo, hi, sign, sumV},
      If[i == j, sumV = 0,
        If[i > j, lo = j; hi = i-1; sign = 1,
           lo = i; hi = j-1; sign = -1];
        sumV = sign * Sum[vdot[p, m]/vRef, {m, lo, hi}];
      ];
      sumV + Vdiff
    ]
  ]
];

Print[""];
Print["=== v_p·(X_i - X_j) / vRef 矩阵 ==="];
edges = {a, b, c, d};
Do[
  Print[""];
  Print["v_", edges[[p]], ":"];
  Do[
    str = "";
    Do[
      If[i != j,
        val = vDotXDifference[edges[[p]], i, j];
        If[val =!= 0,
          str = str <> "  (" <> ToString[i] <> "," <> ToString[j] <> "): " <> ToString[Simplify[val]];
        ];
      ],
      {j, 1, 4}
    ];
    If[str != "", Print[str]],
    {i, 1, 4}
  ],
  {p, 1, 4}
];

(* ========================================================================= *)
(* Step 2: 计算 Δ_ij                                                        *)
(* ========================================================================= *)

(* Δ_ij = -(X_i - X_j)^2 / vRef *)
(* 用 vDotXDifference 的扩展：需要 (X_i - X_j)^2 *)
(* (X_i - X_j)^2 = (V_i-V_j)^2 + 2s_i v_i·(V_i-V_j) - 2s_j v_j·(V_i-V_j) - 2s_i s_j v_i·v_j *)
(* 但 V_i = sum_{k=1}^{i-1} v_k, V_j = sum_{k=1}^{j-1} v_k *)

(* 改用更直接的方法 *)
ClearAll[Delta];
Delta[i_, j_] := Delta[i, j] = Module[
  {vi, vj, si, sj, ViVj2, vDotViVj_i, vDotViVj_j},
  vi = Switch[i, 1,a, 2,b, 3,c, 4,d];
  vj = Switch[j, 1,a, 2,b, 3,c, 4,d];
  si = Switch[i, 1,s1, 2,s2, 3,s3, 4,s4];
  sj = Switch[j, 1,s1, 2,s2, 3,s3, 4,s4];
  If[i == j, Return[0]];
  (* (V_i - V_j)^2 = sum of v_m·v_n for m,n between i and j *)
  (* 用数值方法太复杂，直接用已知公式 *)
  (* 相邻边：|i-j|=1 mod 4 *)
  (* 对边：|i-j|=2 *)
  Module[{diff},
    diff = Mod[Abs[i - j], 4];
    If[diff == 1 || diff == 3,
      (* 相邻边 *)
      Module[{early, late, sEarly, sLate, dot},
        If[(i < j && diff == 1) || (i > j && diff == 3),
          early = i; late = j; sEarly = si; sLate = sj,
          early = j; late = i; sEarly = sj; sLate = si
        ];
        dot = vdot[Switch[early, 1,a, 2,b, 3,c, 4,d], Switch[late, 1,a, 2,b, 3,c, 4,d]]/vRef;
        (1 - sEarly) * sLate * dot
      ],
      (* 对边 *)
      Module[{lo, hi, sLo, sHi, dot},
        If[i < j, lo = i; hi = j; sLo = si; sHi = sj,
           lo = j; hi = i; sLo = sj; sHi = si];
        If[{lo, hi} == {1, 3},
          z*(1 - sLo)*(1 - sHi) + sLo*sHi,
          (1 - sLo)*(1 - sHi) + z*sLo*sHi
        ]
      ]
    ]
  ]
];

Print[""];
Print["=== Δ_ij / vRef ==="];
Do[
  If[i < j,
    val = Delta[i, j];
    Print["Δ_", i, j, " = ", If[val === 0, "0", ToString[Simplify[val]] <> " * vRef"]];
  ],
  {i, 1, 4}, {j, 1, 4}
];

(* ========================================================================= *)
(* Step 3: B_5 和 U                                                          *)
(* ========================================================================= *)

ClearAll[a1, a2, a3, a4, a5];
SA = a1 + a2 + a5;
SB = a3 + a4 + a5;

(* B_5 = Σ_{i<j} α_i α_j c_{ij} Δ_{ij} *)
(* c_{ij} = SB if i,j∈{1,2}; SA if i,j∈{3,4}; α_5 if cross *)
B5 = a1*a2*SB*Delta[1,2] + a3*a4*SA*Delta[3,4] +
     a5*a1*a3*Delta[1,3] + a5*a1*a4*Delta[1,4] +
     a5*a2*a3*Delta[2,3] + a5*a2*a4*Delta[2,4];

U = (a1+a2+a5)*(a3+a4+a5) - a5^2;

Print[""];
Print["=== B_5 / vRef ==="];
Print[Simplify[B5]];

Print[""];
Print["=== U ==="];
Print[U];

(* ========================================================================= *)
(* Step 4: 算子表（28 项）                                                  *)
(* ========================================================================= *)

(* 从 §4 Step 2 的因子化形式构造算子表 *)
(* 每项 = {C_t, v_p_edge, d_p, v_q_edge, d_q} *)
(* d_p = {d1, d2, d3, d4} *)

operatorTerms = {};

(* G1: v_ab*v_cd * (ℓ₂-ℓ₁)·(ℓ₄-ℓ₃) → 4 个 Type I 子项 *)
(* Type I: ℓ_p·ℓ_q → -∂_p·∂_q *)
(* 不适用 (v·D)(v·D) 形式，单独处理 *)
(* G1 展开: vab*vcd*(l2l4 - l2l3 - l1l4 + l1l3) *)
(* O_G1 = -vab*vcd*(∂_2·∂_4 - ∂_2·∂_3 - ∂_1·∂_4 + ∂_1·∂_3) *)

(* G2: v_ab * (C1*B2 + C2*B1) *)
(* C1 = v_d·ℓ₁ - v_d·ℓ₂, B2 = v_c·ℓ₃ + 2v_c·ℓ₄ *)
(* C2 = v_c·ℓ₂ - v_c·ℓ₁, B1 = 2v_d·ℓ₃ + v_d·ℓ₄ *)
(* C1*B2 = (vdl1-vdl2)*(vcl3+2vcl4) *)
(*   = vdl1*vcl3 + 2vdl1*vcl4 - vdl2*vcl3 - 2vdl2*vcl4 *)
(* 4 项: *)
(*   {v_ab, v_d, {+1,0,0,0}, v_c, {0,0,+1,0}}  → v_ab*(v_d·∂₁)(v_c·∂₃) *)
(*   {2v_ab, v_d, {+1,0,0,0}, v_c, {0,0,0,+1}}  → 2v_ab*(v_d·∂₁)(v_c·∂₄) *)
(*   {-v_ab, v_d, {0,+1,0,0}, v_c, {0,0,+1,0}} → -v_ab*(v_d·∂₂)(v_c·∂₃) *)
(*   {-2v_ab, v_d, {0,+1,0,0}, v_c, {0,0,0,+1}} → -2v_ab*(v_d·∂₂)(v_c·∂₄) *)

(* C2*B1 = (vcl2-vcl1)*(2vdl3+vdl4) *)
(*   = 2vcl2*vdl3 + vcl2*vdl4 - 2vcl1*vdl3 - vcl1*vdl4 *)
(* 4 项: *)
(*   {2v_ab, v_c, {0,+1,0,0}, v_d, {0,0,+1,0}} *)
(*   {v_ab, v_c, {0,+1,0,0}, v_d, {0,0,0,+1}} *)
(*   {-2v_ab, v_c, {+1,0,0,0}, v_d, {0,0,+1,0}} *)
(*   {-v_ab, v_c, {+1,0,0,0}, v_d, {0,0,0,+1}} *)

(* G3: v_cd * (D1*A1 + D2*A2) *)
(* D1 = v_a·ℓ₄ - v_a·ℓ₃, A1 = 2v_b·ℓ₁ + v_b·ℓ₂ *)
(* D1*A1 = (val4-val3)*(2vbl1+vbl2) = 2val4*vbl1 + val4*vbl2 - 2val3*vbl1 - val3*vbl2 *)
(* 4 项: *)
(*   {2v_cd, v_a, {0,0,0,+1}, v_b, {+1,0,0,0}} *)
(*   {v_cd, v_a, {0,0,0,+1}, v_b, {0,+1,0,0}} *)
(*   {-2v_cd, v_a, {0,0,+1,0}, v_b, {+1,0,0,0}} *)
(*   {-v_cd, v_a, {0,0,+1,0}, v_b, {0,+1,0,0}} *)

(* D2*A2 = (vbl3-vbl4)*(val1+2val2) = vbl3*val1 + 2vbl3*val2 - vbl4*val1 - 2vbl4*val2 *)
(* 4 项: *)
(*   {v_cd, v_b, {0,0,+1,0}, v_a, {+1,0,0,0}} *)
(*   {2v_cd, v_b, {0,0,+1,0}, v_a, {0,+1,0,0}} *)
(*   {-v_cd, v_b, {0,0,0,+1}, v_a, {+1,0,0,0}} *)
(*   {-2v_cd, v_b, {0,0,0,+1}, v_a, {0,+1,0,0}} *)

(* G4a: v_ac * A1*B1 *)
(* A1*B1 = (2vbl1+vbl2)*(2vdl3+vdl4) = 4vbl1*vdl3 + 2vbl1*vdl4 + 2vbl2*vdl3 + vbl2*vdl4 *)
(*   {4v_ac, v_b, {+1,0,0,0}, v_d, {0,0,+1,0}} *)
(*   {2v_ac, v_b, {+1,0,0,0}, v_d, {0,0,0,+1}} *)
(*   {2v_ac, v_b, {0,+1,0,0}, v_d, {0,0,+1,0}} *)
(*   {v_ac, v_b, {0,+1,0,0}, v_d, {0,0,0,+1}} *)

(* G4b: -v_ad * A1*B2 *)
(* A1*B2 = (2vbl1+vbl2)*(vcl3+2vcl4) = 2vbl1*vcl3 + 4vbl1*vcl4 + vbl2*vcl3 + 2vbl2*vcl4 *)
(*   {-2v_ad, v_b, {+1,0,0,0}, v_c, {0,0,+1,0}} *)
(*   {-4v_ad, v_b, {+1,0,0,0}, v_c, {0,0,0,+1}} *)
(*   {-v_ad, v_b, {0,+1,0,0}, v_c, {0,0,+1,0}} *)
(*   {-2v_ad, v_b, {0,+1,0,0}, v_c, {0,0,0,+1}} *)

(* G4c: -v_bc * A2*B1 *)
(* A2*B1 = (val1+2val2)*(2vdl3+vdl4) = 2val1*vdl3 + val1*vdl4 + 4val2*vdl3 + 2val2*vdl4 *)
(*   {-2v_bc, v_a, {+1,0,0,0}, v_d, {0,0,+1,0}} *)
(*   {-v_bc, v_a, {+1,0,0,0}, v_d, {0,0,0,+1}} *)
(*   {-4v_bc, v_a, {0,+1,0,0}, v_d, {0,0,+1,0}} *)
(*   {-2v_bc, v_a, {0,+1,0,0}, v_d, {0,0,0,+1}} *)

(* G4d: v_bd * A2*B2 *)
(* A2*B2 = (val1+2val2)*(vcl3+2vcl4) = val1*vcl3 + 2val1*vcl4 + 2val2*vcl3 + 4val2*vcl4 *)
(*   {v_bd, v_a, {+1,0,0,0}, v_c, {0,0,+1,0}} *)
(*   {2v_bd, v_a, {+1,0,0,0}, v_c, {0,0,0,+1}} *)
(*   {2v_bd, v_a, {0,+1,0,0}, v_c, {0,0,+1,0}} *)
(*   {4v_bd, v_a, {0,+1,0,0}, v_c, {0,0,0,+1}} *)

(* 总项数: 4(G1) + 8(G2) + 8(G3) + 4+4+4+4(G4) = 36 *)

Print[""];
Print["=== 算子表 ==="];
Print["G1: 4 项 (Type I, ℓ·ℓ)"];
Print["G2: 8 项 (Type II, v_ab 前因子)"];
Print["G3: 8 项 (Type II, v_cd 前因子)"];
Print["G4a: 4 项, G4b: 4 项, G4c: 4 项, G4d: 4 项 (Type III)"];
Print["总计: 36 项"];

(* ========================================================================= *)
(* Step 5: 计算 F_i^(p) 和 H_ij^(pq) [核心]                                *)
(* ========================================================================= *)

(* 系数矩阵 c_ij *)
ClearAll[cMatrix];
Do[
  If[i < j,
    If[i <= 2 && j <= 2, cMatrix[i,j] = SB];
    If[i >= 3 && j >= 3, cMatrix[i,j] = SA];
    If[i <= 2 && j >= 3, cMatrix[i,j] = a5];
  ];
  cMatrix[j,i] = cMatrix[i,j],
  {i, 1, 4}, {j, 1, 4}
];

(* F_i^(p) = -2 α_i Σ_{j≠i} α_j c_{ij} v_p·(X_i - X_j) / vRef *)
(* 返回除以 vRef 后的值 *)
ClearAll[alpha];
alpha[1] = a1; alpha[2] = a2; alpha[3] = a3; alpha[4] = a4;

F[i_, p_] := F[i, p] = -2*alpha[i]*Sum[
  alpha[j] * cMatrix[i,j] * vDotXDifference[p, i, j],
  {j, Complement[{1,2,3,4}, {i}]}
];

(* H_ij^(pq) = 2 α_i α_j c_{ij} (v_p·v_q)/vRef (i≠j) *)
(*           = -2 α_i (v_p·v_q)/vRef Σ_{k≠i} α_k c_{ik} (i=j) *)
H[i_, j_, p_, q_] := H[i, j, p, q] = If[i == j,
  -2*alpha[i]*(vdot[p,q]/vRef)*Sum[alpha[k]*cMatrix[i,k], {k, Complement[{1,2,3,4}, {i}]}],
  2*alpha[i]*alpha[j]*cMatrix[i,j]*(vdot[p,q]/vRef)
];

Print[""];
Print["=== F_i^(p)  (16 个值) ==="];
Do[
  val = Simplify[F[i, edges[[p]]]];
  If[val =!= 0,
    Print["F_", i, "^(v_", edges[[p]], ") = ", val];
  ],
  {i, 1, 4}, {p, 1, 4}
];

(* ========================================================================= *)
(* Step 6: 对每个算子项计算 F^(p), G^(q), H^(pq)                           *)
(* ========================================================================= *)

(* 定义算子项列表 *)
(* 每项: {C_t, v_p_edge, d_vec, v_q_edge, d_vec} *)
(* d_vec = {d1, d2, d3, d4} *)

(* 构建完整的算子项列表 *)
opTerms = {};

(* G2 项 *)
vabVal = vdot[a, b]/vRef;
AppendTo[opTerms, {vabVal, d, {1,0,0,0}, c, {0,0,1,0}}];     (* vdl1*vcl3 *)
AppendTo[opTerms, {2*vabVal, d, {1,0,0,0}, c, {0,0,0,1}}];   (* 2vdl1*vcl4 *)
AppendTo[opTerms, {-vabVal, d, {0,1,0,0}, c, {0,0,1,0}}];    (* -vdl2*vcl3 *)
AppendTo[opTerms, {-2*vabVal, d, {0,1,0,0}, c, {0,0,0,1}}];  (* -2vdl2*vcl4 *)
AppendTo[opTerms, {2*vabVal, c, {0,1,0,0}, d, {0,0,1,0}}];   (* 2vcl2*vdl3 *)
AppendTo[opTerms, {vabVal, c, {0,1,0,0}, d, {0,0,0,1}}];     (* vcl2*vdl4 *)
AppendTo[opTerms, {-2*vabVal, c, {1,0,0,0}, d, {0,0,1,0}}];  (* -2vcl1*vdl3 *)
AppendTo[opTerms, {-vabVal, c, {1,0,0,0}, d, {0,0,0,1}}];    (* -vcl1*vdl4 *)

(* G3 项 *)
vcdVal = vdot[c, d]/vRef;
AppendTo[opTerms, {2*vcdVal, a, {0,0,0,1}, b, {1,0,0,0}}];   (* 2val4*vbl1 *)
AppendTo[opTerms, {vcdVal, a, {0,0,0,1}, b, {0,1,0,0}}];     (* val4*vbl2 *)
AppendTo[opTerms, {-2*vcdVal, a, {0,0,1,0}, b, {1,0,0,0}}];  (* -2val3*vbl1 *)
AppendTo[opTerms, {-vcdVal, a, {0,0,1,0}, b, {0,1,0,0}}];    (* -val3*vbl2 *)
AppendTo[opTerms, {vcdVal, b, {0,0,1,0}, a, {1,0,0,0}}];     (* vbl3*val1 *)
AppendTo[opTerms, {2*vcdVal, b, {0,0,1,0}, a, {0,1,0,0}}];   (* 2vbl3*val2 *)
AppendTo[opTerms, {-vcdVal, b, {0,0,0,1}, a, {1,0,0,0}}];    (* -vbl4*val1 *)
AppendTo[opTerms, {-2*vcdVal, b, {0,0,0,1}, a, {0,1,0,0}}];  (* -2vbl4*val2 *)

(* G4a 项 *)
vacVal = vdot[a, c]/vRef;
AppendTo[opTerms, {4*vacVal, b, {1,0,0,0}, d, {0,0,1,0}}];
AppendTo[opTerms, {2*vacVal, b, {1,0,0,0}, d, {0,0,0,1}}];
AppendTo[opTerms, {2*vacVal, b, {0,1,0,0}, d, {0,0,1,0}}];
AppendTo[opTerms, {vacVal, b, {0,1,0,0}, d, {0,0,0,1}}];

(* G4b 项 *)
vadVal = vdot[a, d]/vRef;
AppendTo[opTerms, {-2*vadVal, b, {1,0,0,0}, c, {0,0,1,0}}];
AppendTo[opTerms, {-4*vadVal, b, {1,0,0,0}, c, {0,0,0,1}}];
AppendTo[opTerms, {-vadVal, b, {0,1,0,0}, c, {0,0,1,0}}];
AppendTo[opTerms, {-2*vadVal, b, {0,1,0,0}, c, {0,0,0,1}}];

(* G4c 项 *)
vbcVal = vdot[b, c]/vRef;
AppendTo[opTerms, {-2*vbcVal, a, {1,0,0,0}, d, {0,0,1,0}}];
AppendTo[opTerms, {-vbcVal, a, {1,0,0,0}, d, {0,0,0,1}}];
AppendTo[opTerms, {-4*vbcVal, a, {0,1,0,0}, d, {0,0,1,0}}];
AppendTo[opTerms, {-2*vbcVal, a, {0,1,0,0}, d, {0,0,0,1}}];

(* G4d 项 *)
vbdVal = vdot[b, d]/vRef;
AppendTo[opTerms, {vbdVal, a, {1,0,0,0}, c, {0,0,1,0}}];
AppendTo[opTerms, {2*vbdVal, a, {1,0,0,0}, c, {0,0,0,1}}];
AppendTo[opTerms, {2*vbdVal, a, {0,1,0,0}, c, {0,0,1,0}}];
AppendTo[opTerms, {4*vbdVal, a, {0,1,0,0}, c, {0,0,0,1}}];

Print[""];
Print["=== 算子项数 ==="];
Print[Length[opTerms], " 项 (Type II+III，不含 Type I 的 4 项)"];

(* 对每项计算 F^(p), G^(q), H^(pq) *)
ClearAll[FG, HHterm];
Do[
  {Ct, vp, dp, vq, dq} = opTerms[[t]];
  (* F^(p) = Σ d_i F_i^(vp) *)
  FGp = Sum[dp[[i]] * F[i, vp], {i, 1, 4}];
  (* G^(q) = Σ d_j F_j^(vq) *)
  FGq = Sum[dq[[j]] * F[j, vq], {j, 1, 4}];
  (* H^(pq) = Σ d_i d_j H_ij^(vp, vq) *)
  HHpq = Sum[dp[[i]] * dq[[j]] * H[i, j, vp, vq], {i, 1, 4}, {j, 1, 4}];
  FG[t] = Simplify[Ct * FGp * FGq];
  HHterm[t] = Simplify[Ct * HHpq];
  ,
  {t, 1, Length[opTerms]}
];

(* Type I (G1) 的贡献 *)
(* G1 = v_ab*v_cd*(ℓ₂-ℓ₁)·(ℓ₄-ℓ₃) *)
(* ℓ→i∂_X: O_G1 = -v_ab*v_cd*(∂_2-∂_1)·(∂_4-∂_3) *)
(* Type I 算子作用于 B_5^n: *)
(* (∂_i·∂_j) B_5^n = n(n-1) B_5^{n-2} (∂_i B_5)·(∂_j B_5) + n B_5^{n-1} (∂_i·∂_j B_5) *)

(* 对于 Type I，有 4 个 (i,j) 对 *)
typeI = {{2,4,1}, {2,3,-1}, {1,4,-1}, {1,3,1}};  (* {i, j, sign} *)

(* F_i^(I) = ∂_i B_5 = Σ_{k≠i} 2α_i α_k c_{ik} (X_i - X_k) *)
(* 注意：Δ_ij = -(X_i-X_j)^2，所以 B_5 = Σ α_i α_j c_{ij} Δ_ij = -Σ α_i α_j c_{ij} (X_i-X_j)^2 *)
(* ∂_i B_5 = -2α_i Σ_{k≠i} α_k c_{ik} (X_i-X_k) *)
(* 恰好等于 F_i^(p) 去掉 v_p·! 即 F_i^(I) = ∂_i B_5 = -2α_i Σ α_k c_{ik} (X_i-X_k) [向量] *)

(* 但 F_i^(p) = v_p·∂_i B_5，所以 ∂_i B_5 = Σ_{p} F_i^(vp) * (某个基向量) *)
(* 对于 Type I，我们需要 F_i^(I) · F_j^(I) = (∂_i B_5)·(∂_j B_5) *)
(* 和 H_ij^(I) = ∂_i·∂_j B_5 = D * 2α_i α_j c_{ij} (i≠j) *)

(* 简化：直接计算 *)
(* (∂_i B_5)·(∂_j B_5) = 4α_i α_j Σ_{k≠i} Σ_{l≠j} α_k α_l c_{ik} c_{jl} (X_i-X_k)·(X_j-X_l) *)
(* (X_i-X_k)·(X_j-X_l) = 1/2 [Δ_{kl} + Δ_{ij} - Δ_{ik} - Δ_{jl}] *)
(* 这非常复杂，交给 montopar 代码处理 *)

(* 对于 Type I 的 Pnum 贡献： *)
(* Pnum_I^(1) = (2-3ε) Σ C_t F_I^(t) G_I^(t) *)
(* Pnum_I^(2) = -Σ C_t H_I^(t) *)
(* 其中 C_t = ±v_ab*v_cd *)

Print[""];
Print["=== Pnum 贡献汇总 ==="];

(* 计算 Pnum^(1) 和 Pnum^(2) 来自 Type II+III *)
Pnum1_IIIII = Simplify[Sum[FG[t], {t, 1, Length[opTerms]}]];
Pnum2_IIIII = Simplify[-Sum[HHterm[t], {t, 1, Length[opTerms]}]];

Print["Pnum^(1) [Type II+III, 不含 (2-3ε) 因子]: "];
Print[Pnum1_IIIII];

Print[""];
Print["Pnum^(2) [Type II+III]: "];
Print[Pnum2_IIIII];

(* ========================================================================= *)
(* Step 7: C 因子                                                            *)
(* ========================================================================= *)

Print[""];
Print["=============================================="];
Print["C 和 J 汇总"];
Print["=============================================="];

Print[""];
Print["=== C 因子 ==="];
Print["C = ① color × (4π α_s)^3 × Γ(1-3ε)/(1024 π^{6-3ε}) × (1-3ε) × (-2v_{14}+i0)^{3ε}"];
Print[""];
Print["① color = SUNTr[T^{A1}T^{A2}T^{A3}T^{A4}] f^{A1 A2 B} f^{A3 A4 B}"];
Print["  (需 FeynCalc 求值，轮廓排序: A1, A2, A3, A4)"];

(* ========================================================================= *)
(* Step 8: J 积分                                                            *)
(* ========================================================================= *)

Print[""];
Print["=== J 积分 ==="];
Print["J = ∫₀¹ d⁴s ∫ d⁵α δ(1-Σα) (α₁α₂α₃α₄α₅)^{-ε} / U^{1+2ε}"];
Print["    × [Pnum^(1)_total / B_5^{3-3ε} + Pnum^(2)_total / B_5^{2-3ε}]"];

Print[""];
Print["其中:"];
Print["  U = (α₁+α₂+α₅)(α₃+α₄+α₅) - α₅²"];
Print["  B_5 = ", Simplify[B5]];
Print["  Pnum^(1)_total = (2-3ε) × [Pnum^(1)_IIIII + Pnum^(1)_I]"];
Print["  Pnum^(2)_total = Pnum^(2)_IIIII + Pnum^(2)_I"];

Print[""];
Print["  Pnum^(1)_IIIII (不含 (2-3ε) 因子) = 见上方"];
Print["  Pnum^(2)_IIIII = 见上方"];
Print["  Pnum^(1)_I, Pnum^(2)_I = Type I (G1) 的贡献，需额外计算"];

Print[""];
Print["v07 完成"];