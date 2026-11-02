(* v02-color-factors.m *)
(* ========================================================================= *)
(* 验证：各边分布的颜色因子                                              *)
(* 方法：FeynCalc SUNSimplify，使用与 montopar 模板相同的约定              *)
(*                                                                           *)
(* 约定：                                                                    *)
(*   SUNT[a1, a2, a3, a4] = Tr[T^{a1} T^{a2} T^{a3} T^{a4}]                  *)
(*   索引 a1,a2,a3,a4 按轮廓顺序排列（边 1→2→3→4，同边内按 s 参数）       *)
(*   TGV A: SUNF[a,b,B]，TGV B: SUNF[c,d,B]                                  *)
(*   配对固定：TGV A 腿 (a,b)，TGV B 腿 (c,d)                               *)
(* ========================================================================= *)

Print["=============================================="];
Print["v02: 颜色因子验证"];
Print["=============================================="];

<< FeynCalc`

(*
  N=4 四边形，腿仅在边 3,4 上。

  配对固定：
    腿 1,2 → TGV A（边 a,b）
    腿 3,4 → TGV B（边 c,d）

  三种配对：
    S 型 (33)(44): a=3,b=3,c=4,d=4 → 腿 1,2 在边 3，腿 3,4 在边 4
    T 型 (34)(34): a=3,b=4,c=3,d=4 → 腿 1,3 在边 3，腿 2,4 在边 4
    U 型 (34)(43): a=3,b=4,c=4,d=3 → 腿 1,4 在边 3，腿 2,3 在边 4
*)

Print[""];
Print["=== 边 3,4 上的分布 ==="];

(* S 型 (3,3,4,4)：轮廓排序 a1,a2,a3,a4 *)
(* SUNF[A1,A2,B] * SUNF[A3,A4,B] * SUNT[A1,A2,A3,A4] *)
(* 注意：索引 A1,A2,A3,A4 同时出现在 SUNT 和 SUNF 中，SUNSimplify 会缩并 *)

colorS = SUNSimplify[
  SUNF[A1, A2, B] * SUNF[A3, A4, B] * SUNT[A1, A2, A3, A4],
  SUNNToCACF -> True
];
Print["S 型 (3,3,4,4): ", colorS];

(* T 型 (3,4,3,4)：轮廓 a1 在边 3，a3 在边 3，a2 在边 4，a4 在边 4 *)
(* 排序：A1, A3, A2, A4 *)
(* SUNF 仍用 A1,A2 和 A3,A4（TGV 配对不随轮廓排序改变）*)
colorT = SUNSimplify[
  SUNF[A1, A2, B] * SUNF[A3, A4, B] * SUNT[A1, A3, A2, A4],
  SUNNToCACF -> True
];
Print["T 型 (3,4,3,4): ", colorT];

(* U 型 (3,4,4,3)：边 3: A1, A4；边 4: A2, A3 *)
(* 排序：A1, A4, A2, A3 *)
colorU = SUNSimplify[
  SUNF[A1, A2, B] * SUNF[A3, A4, B] * SUNT[A1, A4, A2, A3],
  SUNNToCACF -> True
];
Print["U 型 (3,4,4,3): ", colorU];

(* ========================================================================= *)
(* 与文档对照                                                               *)
(* ========================================================================= *)

Print[""];
Print["=== 对照 general-derivation.md §3b ==="];

check[name_, fc_, expected_] := Module[{diff},
  diff = SUNSimplify[fc - expected, SUNNToCACF -> True];
  If[diff === 0,
    Print["  ", name, ": \[Checkmark]"];
    ,
    Print["  ", name, ": \[WarningSign] (差异: ", diff, ")"]
  ]
];

check["S 型 (3,3,4,4) = -CA^2*CF/4", colorS, -CA^2*CF/4];
check["T 型 (3,4,3,4) = +CA^2*CF/4", colorT, CA^2*CF/4];
check["U 型 (3,4,4,3) = 0", colorU, 0];

(* ========================================================================= *)
(* Jacobi 恒等式                                                            *)
(* ========================================================================= *)

Print[""];
Print["=== Jacobi 恒等式: color_S + color_T + color_U = 0 ==="];
jacobi = SUNSimplify[colorS + colorT + colorU, SUNNToCACF -> True];
Print["  ", If[jacobi === 0, "\[Checkmark]", "\[WarningSign] 结果: " <> ToString[jacobi]]];

(* ========================================================================= *)
(* 扩展到所有边分布                                                         *)
(* ========================================================================= *)

Print[""];
Print["=== 其他边分布（供参考）==="];

(* 辅助函数 *)
computeColor[contourOrder_List] := SUNSimplify[
  SUNF[A1, A2, B] * SUNF[A3, A4, B] * SUNT @@ contourOrder,
  SUNNToCACF -> True
];

(* 所有边分布及其轮廓排序 *)
distributions = {
  {"(1,1,2,2)", {A1, A2, A3, A4}},
  {"(1,2,1,2)", {A1, A3, A2, A4}},
  {"(1,2,2,1)", {A1, A4, A2, A3}},
  {"(1,1,3,3)", {A1, A2, A3, A4}},
  {"(1,3,1,3)", {A1, A3, A2, A4}},
  {"(1,3,3,1)", {A1, A4, A2, A3}},
  {"(2,2,3,3)", {A1, A2, A3, A4}},
  {"(2,3,2,3)", {A1, A3, A2, A4}},
  {"(2,3,3,2)", {A1, A4, A2, A3}},
  {"(1,2,3,4)", {A1, A2, A3, A4}},
  {"(1,3,2,4)", {A1, A3, A2, A4}},
  {"(1,4,2,3)", {A1, A4, A2, A3}}
};

TableForm[
  Table[
    {name, computeColor[order]},
    {name, order} = dist
  ],
  TableHeadings -> {None, {"边分布", "颜色因子"}}
];

Print[""];
Print["v02 完成"];