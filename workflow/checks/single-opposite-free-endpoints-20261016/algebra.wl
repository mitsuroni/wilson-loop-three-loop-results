Module[{zz, ll, aa, bb, ff, nu, th, op, opg, vals, poly, mm, jj},
 th[x_] := zz D[x, zz];
 op[x_] := th[th[x]] - ll th[x] + zz (th[th[x]] + (2 + aa + bb) th[x] + (1 + aa) (1 + bb) x);
 opg[x_] := th[th[x]] - ll th[x] + zz (th[th[x]] - (aa + bb + 2 ll) th[x] + (aa + ll) (bb + ll) x);
 nu = 1 + aa + bb + ll;
 vals = {FullSimplify[PowerExpand[(1 + zz)^nu op[(1 + zz)^(-nu) ff[zz]]] - opg[ff[zz]]] === 0};
 Do[
  poly = Sum[Binomial[mm - 1, jj] zz^jj Factorial[mm - jj - 1]^2/
    (Pochhammer[1 + aa + jj, mm - jj] Pochhammer[1 + bb + jj, mm - jj]), {jj, 0, mm - 1}];
  AppendTo[vals, Together[(op[poly] /. ll -> mm) - zz^mm] === 0];
  AppendTo[vals, Together[(op[poly] /. ll -> mm) + zz^mm] =!= 0],
  {mm, 1, 4}];
 AppendTo[vals, FullSimplify[Gamma[-ll]/Gamma[1 - ll] + 1/ll] === 0];
 AppendTo[vals, FullSimplify[PolyGamma[0, -ll] - PolyGamma[0, 1 - ll] - 1/ll] === 0];
 <|"Passed" -> And @@ vals, "Checks" -> vals, "Count" -> Length[vals], "Formal" -> "not_run"|>
]
