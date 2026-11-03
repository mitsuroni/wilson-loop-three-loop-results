Module[{eps, ell, center, width, widths, d1, d2, d3, pref, interval, product, variance, want4, checks},
 widths = {d1, d2, d3};
 pref = Normal@Series[Exp[-3 EulerGamma eps] Gamma[1-eps]^3 Sech[3 eps ell/2], {eps,0,4}];
 interval = Normal@Series[Exp[eps center] Sinh[eps width/2]/(eps width/2), {eps,0,4}];
 product = Normal@Series[Exp[eps center] Times@@(Sinh[eps #/2]/(eps #/2)& /@ widths), {eps,0,4}];
 variance = Total[widths^2]/24;
 want4 = center^4/24 + center^2 variance/2 + variance^2/2 - Total[widths^4]/2880;
 checks = <|
  "prefactor_eps3" -> (FullSimplify[Coefficient[pref,eps,3]-Zeta[3]]===0),
  "prefactor_eps4" -> (FullSimplify[Coefficient[pref,eps,4]-(3 Zeta[4]/4+9 Zeta[2]^2/8-27 Zeta[2] ell^2/16+135 ell^4/128)]===0),
  "interval_eps3" -> (Expand[Coefficient[interval,eps,3]-(center^3/6+center width^2/24)]===0),
  "interval_eps4" -> (Expand[Coefficient[interval,eps,4]-(center^4/24+center^2 width^2/48+width^4/1920)]===0),
  "three_intervals_eps3" -> (Expand[Coefficient[product,eps,3]-(center^3/6+center variance)]===0),
  "three_intervals_eps4" -> (Expand[Coefficient[product,eps,4]-want4]===0),
  "missing_interval_moment_rejected" -> (Expand[Coefficient[interval,eps,4]-center^4/24]=!=0)
 |>;
 Append[checks,"all_pass" -> And@@Values[checks]]
]
