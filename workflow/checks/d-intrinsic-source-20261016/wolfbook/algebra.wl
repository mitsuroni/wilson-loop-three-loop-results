Module[{t, s, z, x, y, ee, ll, vv, hh, qq, rr, ss, aa, bb,
  eta, kap, u, v, ell, assum, cross, cfun, powers, logs, rect,
  jac, norm, positive, negative, checks},
  assum = 0 < t < 1 && 0 < s < 1 && z > 0 && 0 < x < 1 && 0 < y < 1;
  hh = (1 - t) (1 - s) + z t s;
  qq = hh/Sqrt[z]; rr = t (1 - s); ss = (1 - t) s;
  aa = Sqrt[z] t; bb = (1 - s)/Sqrt[z];
  eta = Log[(1 - t)/(z s)]/2;
  kap = (1 + z) Sqrt[ss/z];
  v = Log[qq/Sqrt[rr]]; ell = Log[z t/(1 - s)]/2;
  u = Log[t (1 - t) s (1 - s)]/2 + Log[qq];
  cross = {z -> 1/z, t -> 1 - s, s -> 1 - t};
  cfun[f_] := f /. cross;
  powers = FullSimplify[
    (Exp[ee ll] - Exp[ee vv]) (Exp[ee vv] - Exp[-ee ll]) -
    2 Exp[ee vv] (Cosh[ee ll] - Cosh[ee vv]),
    Element[{ee, ll, vv}, Reals]];
  logs = FullSimplify[PowerExpand[Log[rr] + v + Log[ss]/2 - u], assum];
  rect = FullSimplify[
    (((1 - (1 - (1 - t) x)) (1 - s y) + z (1 - (1 - t) x) s y)/Sqrt[z]) -
    Sqrt[ss] (Exp[eta] x + Exp[-eta] y - kap x y), assum];
  jac = FullSimplify[-Det[{{D[1 - (1 - t) x, x], D[1 - (1 - t) x, y]},
      {D[s y, x], D[s y, y]}}] - ss, assum];
  norm = FullSimplify[(1 + z)/Sqrt[z] Sqrt[ss] - kap, assum];
  positive = {
    "generic power difference" -> (powers === 0),
    "endpoint product" -> (FullSimplify[aa bb - rr, assum] === 0),
    "upper endpoint logarithm" -> (FullSimplify[PowerExpand[Log[aa] - Log[rr]/2 - ell], assum] === 0),
    "lower endpoint logarithm" -> (FullSimplify[PowerExpand[Log[bb] - Log[rr]/2 + ell], assum] === 0),
    "middle endpoint logarithm" -> (FullSimplify[PowerExpand[Log[qq] - Log[rr]/2 - v], assum] === 0),
    "complete exponent" -> (logs === 0),
    "rectangle coordinate" -> (rect === 0),
    "absolute rectangle Jacobian" -> (jac === 0),
    "rectangle normalization" -> (norm === 0),
    "crossing q" -> (FullSimplify[cfun[qq] - qq, assum] === 0),
    "crossing R" -> (FullSimplify[cfun[rr] - rr, assum] === 0),
    "crossing S" -> (FullSimplify[cfun[ss] - ss, assum] === 0),
    "crossing u" -> (FullSimplify[PowerExpand[cfun[u] - u], assum] === 0),
    "crossing v" -> (FullSimplify[PowerExpand[cfun[v] - v], assum] === 0),
    "crossing ell" -> (FullSimplify[PowerExpand[cfun[ell] + ell], assum] === 0),
    "crossing eta" -> (FullSimplify[PowerExpand[cfun[eta] + eta], assum] === 0),
    "crossing kappa" -> (FullSimplify[cfun[kap] - kap, assum] === 0),
    "crossing Cauchy weight" -> (FullSimplify[
      cfun[(1 + z)^2/(((1 + z) t - 1) ((1 + z) s - 1))] -
      (1 + z)^2/(((1 + z) t - 1) ((1 + z) s - 1)), assum] === 0),
    "crossing rectangle exchange" -> (FullSimplify[
      Exp[-eta] y + Exp[eta] x - kap y x -
      (Exp[eta] x + Exp[-eta] y - kap x y), assum] === 0)
  };
  negative = {
    "missing factor two rejected" -> TrueQ[
      (Exp[ee ll] - Exp[ee vv]) (Exp[ee vv] - Exp[-ee ll]) !=
      Exp[ee vv] (Cosh[ee ll] - Cosh[ee vv]) /.
      {ee -> 1, ll -> Log[2], vv -> Log[3]}],
    "wrong bilinear sign rejected" -> TrueQ[
      (Sqrt[ss] 2 kap x y /. {t -> 1/3, s -> 2/5, z -> 7/5, x -> 1/2, y -> 1/3}) != 0],
    "missing rectangle Jacobian rejected" -> TrueQ[(ss /. {t -> 1/3, s -> 2/5}) != 1],
    "unflipped eta rejected" -> TrueQ[
      (Exp[eta] x + Exp[-eta] y - Exp[-eta] x - Exp[eta] y /.
      {t -> 1/3, s -> 2/5, z -> 7/5, x -> 1/2, y -> 1/3}) != 0]
  };
  checks = Join[positive, negative];
  Print[Grid[Prepend[({First[#], Last[#]} & /@ checks), {"Check", "Passed"}], Frame -> All]];
  And @@ (Last /@ checks)
]
