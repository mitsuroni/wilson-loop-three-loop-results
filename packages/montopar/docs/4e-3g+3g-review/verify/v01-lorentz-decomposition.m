(* v01-lorentz-decomposition.m *)
(* ========================================================================= *)
(* 小任务：V_{2TGV} 有效顶点的 Lorentz 分解                              *)
(* 从三胶子顶点出发，独立推导，用符号替换做系数提取，与文档对照           *)
(* ========================================================================= *)

Print["=============================================="];
Print["v01: V_{2TGV} Lorentz 分解"];
Print["=============================================="];

<< FeynCalc`

(* ========================================================================= *)
(* Step 1: FeynCalc 直接计算                                                *)
(* ========================================================================= *)

ClearAll[Gamma3g];
Gamma3g[mu_, nu_, rho_, k1_, k2_, k3_] :=
  MTD[mu, nu] * FVD[k1 - k2, rho] +
  MTD[nu, rho] * FVD[k2 - k3, mu] +
  MTD[rho, mu] * FVD[k3 - k1, nu];

tgvA = Contract[FVD[a, mu1] * FVD[b, mu2] * Gamma3g[mu1, mu2, nu5, -l1, -l2, l1 + l2]];
tgvB = Contract[FVD[c, mu3] * FVD[d, mu4] * Gamma3g[mu3, mu4, nu6, -l3, -l4, l3 + l4]];
fullLorentz = Contract[MTD[nu5, nu6] * tgvA * tgvB];
lorentzFC = ExpandScalarProduct[fullLorentz];

(* 转换为外部可读形式 *)
lorentzExt = FCE[lorentzFC];
Print["FeynCalc 直接结果: ", lorentzExt];

(* ========================================================================= *)
(* Step 2: 用 SPD 构造 9 项并验证                                           *)
(* ========================================================================= *)

(* 9 项分解 *)
A1B1 = SPD[a, b] * SPD[c, d] * (SPD[l2, l4] - SPD[l2, l3] - SPD[l1, l4] + SPD[l1, l3]);
A1B2 = SPD[a, b] * (SPD[d, l2] - SPD[d, l1]) * (-SPD[c, l3] - 2*SPD[c, l4]);
A1B3 = SPD[a, b] * (SPD[c, l2] - SPD[c, l1]) * (2*SPD[d, l3] + SPD[d, l4]);
A2B1 = SPD[c, d] * (SPD[b, l4] - SPD[b, l3]) * (-SPD[a, l1] - 2*SPD[a, l2]);
A2B2 = SPD[b, d] * (-SPD[a, l1] - 2*SPD[a, l2]) * (-SPD[c, l3] - 2*SPD[c, l4]);
A2B3 = SPD[b, c] * (-SPD[a, l1] - 2*SPD[a, l2]) * (2*SPD[d, l3] + SPD[d, l4]);
A3B1 = SPD[c, d] * (SPD[a, l4] - SPD[a, l3]) * (2*SPD[b, l1] + SPD[b, l2]);
A3B2 = SPD[a, d] * (2*SPD[b, l1] + SPD[b, l2]) * (-SPD[c, l3] - 2*SPD[c, l4]);
A3B3 = SPD[a, c] * (2*SPD[b, l1] + SPD[b, l2]) * (2*SPD[d, l3] + SPD[d, l4]);

manual9 = A1B1 + A1B2 + A1B3 + A2B1 + A2B2 + A2B3 + A3B1 + A3B2 + A3B3;
manual9Ext = ExpandScalarProduct[manual9] // FCE;

(* 比较：将两个表达式都转为 FCI 内部形式 *)
diff9 = ExpandScalarProduct[fullLorentz - (manual9 // FCI)];
Print[""];
Print["9项手动 vs FeynCalc: ", If[diff9 === 0, "✓ 一致", "✗ 不一致: " <> ToString[FCE[diff9]]]];

(* ========================================================================= *)
(* Step 3: 用符号替换提取系数                                              *)
(* ========================================================================= *)

(* 将 SPD[a,b] 等替换为可做 Coefficient 的符号 *)
toSyms = {
  SPD[a, b] -> vab, SPD[a, c] -> vac, SPD[a, d] -> vad,
  SPD[b, c] -> vbc, SPD[b, d] -> vbd, SPD[c, d] -> vcd,
  SPD[a, l1] -> val1, SPD[a, l2] -> val2, SPD[a, l3] -> val3, SPD[a, l4] -> val4,
  SPD[b, l1] -> vbl1, SPD[b, l2] -> vbl2, SPD[b, l3] -> vbl3, SPD[b, l4] -> vbl4,
  SPD[c, l1] -> vcl1, SPD[c, l2] -> vcl2, SPD[c, l3] -> vcl3, SPD[c, l4] -> vcl4,
  SPD[d, l1] -> vdl1, SPD[d, l2] -> vdl2, SPD[d, l3] -> vdl3, SPD[d, l4] -> vdl4,
  SPD[l1, l1] -> l1l1, SPD[l2, l2] -> l2l2, SPD[l3, l3] -> l3l3, SPD[l4, l4] -> l4l4,
  SPD[l1, l2] -> l1l2, SPD[l1, l3] -> l1l3, SPD[l1, l4] -> l1l4,
  SPD[l2, l3] -> l2l3, SPD[l2, l4] -> l2l4, SPD[l3, l4] -> l3l4
};

symFC = lorentzExt /. toSyms;
Print[""];
Print["符号化后项数: ", Length[symFC]];

(* 提取各 v_{ij} 前因子的系数 *)
coeffVabVcd = Coefficient[symFC, vab*vcd];
coeffVab = Coefficient[symFC - coeffVabVcd*vab*vcd, vab];  (* 不含 vcd 的 vab 项 *)
coeffVcd = Coefficient[symFC - coeffVabVcd*vab*vcd - coeffVab*vab, vcd];  (* 不含 vab 的 vcd 项 *)
coeffVac = Coefficient[symFC, vac];
coeffVad = Coefficient[symFC, vad];
coeffVbc = Coefficient[symFC, vbc];
coeffVbd = Coefficient[symFC, vbd];

(* 验证分组重构 *)
reconstructed = coeffVabVcd*vab*vcd + coeffVab*vab + coeffVcd*vcd +
  coeffVac*vac + coeffVad*vad + coeffVbc*vbc + coeffVbd*vbd;

diffRecon = Expand[symFC - reconstructed];
Print["分组重构一致: ", If[diffRecon === 0, "✓", "✗ 剩余: " <> ToString[diffRecon]]];

(* ========================================================================= *)
(* Step 4: 展示 FeynCalc 直接结果的正确分组                                 *)
(* ========================================================================= *)

Print[""];
Print["=============================================="];
Print["FeynCalc 直接结果的正确分组"];
Print["=============================================="];

Print[""];
Print["G1 (v_ab * v_cd 前因子):"];
Print["  ", coeffVabVcd];

Print[""];
Print["G2 (v_ab 前因子，不含 v_cd):"];
Print["  ", coeffVab];

Print[""];
Print["G3 (v_cd 前因子，不含 v_ab):"];
Print["  ", coeffVcd];

Print[""];
Print["G4a (v_ac 前因子):"];
Print["  ", coeffVac];

Print[""];
Print["G4b (v_ad 前因子):"];
Print["  ", coeffVad];

Print[""];
Print["G4c (v_bc 前因子):"];
Print["  ", coeffVbc];

Print[""];
Print["G4d (v_bd 前因子):"];
Print["  ", coeffVbd];

(* ========================================================================= *)
(* Step 5: 与 effective-vertices.md §3b 对照（用符号形式）                  *)
(* ========================================================================= *)

Print[""];
Print["=============================================="];
Print["对照 effective-vertices.md §3b"];
Print["=============================================="];

ClearAll[vab, vac, vad, vbc, vbd, vcd];
ClearAll[val1, val2, val3, val4, vbl1, vbl2, vbl3, vbl4];
ClearAll[vcl1, vcl2, vcl3, vcl4, vdl1, vdl2, vdl3, vdl4];
ClearAll[l1l1, l2l2, l3l3, l4l4, l1l2, l1l3, l1l4, l2l3, l2l4, l3l4];

(* effective-vertices.md 声称的形式 *)
G1doc = vab * vcd * (l2l4 - l2l3 - l1l4 + l1l3);

G2doc = vab * (
  2*vcl1*vdl2 - 2*vcl2*vdl1 +
  (vcl2 - vcl1)*vdl3 +
  (vdl1 - vdl2)*vcl4
);

G3doc = 2*vcd * (
  (-val1 - 2*val2)*(vbl4 - vbl3) +
  (2*vbl1 + vbl2)*(val4 - val3)
);

G4aDoc = vac * (2*vbl1 + vbl2) * (vdl3 - vdl1 - vdl2);
G4bDoc = vad * (2*vbl1 + vbl2) * (vcl1 + vcl2 - vcl4);
G4cDoc = -vbc * (val1 + 2*val2) * (vdl3 - vdl1 - vdl2);
G4dDoc = -vbd * (val1 + 2*val2) * (vcl1 + vcl2 - vcl4);

(* 对照 *)
check[fc_, doc_, label_] := Module[{diff},
  diff = Expand[fc - doc];
  If[diff === 0,
    Print["  ", label, ": ✓"],
    Print["  ", label, ": ✗"];
    Print["    差异: ", diff]
  ]
];

check[coeffVabVcd, Expand[G1doc / (vab*vcd)] * vab * vcd, "G1 (v_ab*v_cd)"];

(* G1 的系数应该是 l2l4 - l2l3 - l1l4 + l1l3 *)
(* 文档的 G1 = vab*vcd*(l2l4-l2l3-l1l4+l1l3) *)
(* 所以 coeffVabVcd 应该等于 l2l4 - l2l3 - l1l4 + l1l3 *)
check[coeffVabVcd, l2l4 - l2l3 - l1l4 + l1l3, "G1 系数"];

(* G2 *)
check[coeffVab, Expand[G2doc / vab], "G2 (v_ab)"];

(* G3 *)
check[coeffVcd, Expand[G3doc / vcd], "G3 (v_cd)"];

(* G4 *)
check[coeffVac, Expand[G4aDoc / vac], "G4a (v_ac)"];
check[coeffVad, Expand[G4bDoc / vad], "G4b (v_ad)"];
check[coeffVbc, Expand[G4cDoc / vbc], "G4c (v_bc)"];
check[coeffVbd, Expand[G4dDoc / vbd], "G4d (v_bd)"];

(* ========================================================================= *)
(* Step 6: 数值交叉验证                                                      *)
(* ========================================================================= *)

Print[""];
Print["=============================================="];
Print["数值交叉验证（5 组随机运动学）"];
Print["=============================================="];

SeedRandom[42];
randomVec[] := RandomReal[{-1, 1}, 4];
numDot[v1_, v2_] := Sum[v1[[i]]*v2[[i]], {i, 1, 4}];

allPass = True;
Do[
  rv = {randomVec[], randomVec[], randomVec[], randomVec[]};
  rp = {randomVec[], randomVec[], randomVec[], randomVec[]};

  numRules = {
    vab -> numDot[rv[[1]], rv[[2]]], vac -> numDot[rv[[1]], rv[[3]]],
    vad -> numDot[rv[[1]], rv[[4]]], vbc -> numDot[rv[[2]], rv[[3]]],
    vbd -> numDot[rv[[2]], rv[[4]]], vcd -> numDot[rv[[3]], rv[[4]]],
    val1 -> numDot[rv[[1]], rp[[1]]], val2 -> numDot[rv[[1]], rp[[2]]],
    val3 -> numDot[rv[[1]], rp[[3]]], val4 -> numDot[rv[[1]], rp[[4]]],
    vbl1 -> numDot[rv[[2]], rp[[1]]], vbl2 -> numDot[rv[[2]], rp[[2]]],
    vbl3 -> numDot[rv[[2]], rp[[3]]], vbl4 -> numDot[rv[[2]], rp[[4]]],
    vcl1 -> numDot[rv[[3]], rp[[1]]], vcl2 -> numDot[rv[[3]], rp[[2]]],
    vcl3 -> numDot[rv[[3]], rp[[3]]], vcl4 -> numDot[rv[[3]], rp[[4]]],
    vdl1 -> numDot[rv[[4]], rp[[1]]], vdl2 -> numDot[rv[[4]], rp[[2]]],
    vdl3 -> numDot[rv[[4]], rp[[3]]], vdl4 -> numDot[rv[[4]], rp[[4]]],
    l1l1 -> numDot[rp[[1]], rp[[1]]], l2l2 -> numDot[rp[[2]], rp[[2]]],
    l3l3 -> numDot[rp[[3]], rp[[3]]], l4l4 -> numDot[rp[[4]], rp[[4]]],
    l1l2 -> numDot[rp[[1]], rp[[2]]], l1l3 -> numDot[rp[[1]], rp[[3]]],
    l1l4 -> numDot[rp[[1]], rp[[4]]], l2l3 -> numDot[rp[[2]], rp[[3]]],
    l2l4 -> numDot[rp[[2]], rp[[4]]], l3l4 -> numDot[rp[[3]], rp[[4]]]
  };

  nFC = symFC /. numRules;
  nDoc = (G1doc + G2doc + G3doc + G4aDoc + G4bDoc + G4cDoc + G4dDoc) /. numRules;
  relDiff = Abs[nFC - nDoc] / Max[Abs[nFC], 10^-15];

  status = If[relDiff < 10^-12, "PASS", "FAIL"];
  If[relDiff >= 10^-12, allPass = False];
  Print["  试验 ", k, ": ", status, " (diff = ", ScientificForm[relDiff], ")"],
  {k, 1, 5}
];

(* ========================================================================= *)
(* Step 7: 结论                                                              *)
(* ========================================================================= *)

Print[""];
Print["=============================================="];
Print["结论"];
Print["=============================================="];

Print[""];
If[allPass && diff9 === 0,
  Print["✓ 全部通过：V_{2TGV} Lorentz 分解与 effective-vertices.md §3b 一致"],
  Print["✗ 存在差异，effective-vertices.md §3b 需要修正。"];
  Print["  正确形式见上方 Step 4 的 FeynCalc 输出。"]
];

Print[""];
Print["v01 完成"];