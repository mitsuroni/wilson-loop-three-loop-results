# Three-loop Wilson-loop exponent: four-color analytic results

The three-loop exponent is resolved into four color channels. This article gives their complete Laurent coefficients through the finite term, the renormalized finite functions, and the exact connection between the scale-stripped sum and the physical result. The coefficients are written in the established $M,N$ functions and zeta values.

**Read:** [finite functions](#finite-functions) · [all-laurent coefficients](#all-laurent-coefficients) · [function definitions](#function-definitions) · [exact data](../../benchmarks/current/results/3loop.json) · [reproduce](../../docs/reproduction.md).

## Quantity and normalization

Set $D=4-2\epsilon$, $\lambda=\alpha_s^{\overline{\mathrm{MS}}}/\pi$ and $T=T_Fn_f$, with $T_F=1/2$. The renormalized connected exponent is

$$\log W_R=\sum_{n\ge1}\lambda^nF_n,\qquad
F_3=C_F\left[C_A^2F_{3,AA}+C_ATF_{3,Af}+C_FTF_{3,Ff}+T^2F_{3,ff}\right].$$

The same color decomposition applies to $U_3$, $B_3$ and $P_3$. An entry labelled $AA$, for example, is the coefficient after removing $C_FC_A^2$. These are coefficients of the connected logarithm.

For the null quadrilateral, set $s=2v_1\!\cdot v_2$, $t=2v_1\!\cdot v_4$ and $z=s/t$. Use $\rho=\log z$ and $L=(L_s+L_t)/2$, so $L_s=L+\rho/2$ and $L_t=L-\rho/2$. The original source domain, signs and Euler factors follow the [single physical convention](../../workflow/current-conventions.md). The displayed real-function representation uses $z>0$.

The calculation restores the external scale, adds coupling counterterms and separates the pole and finite parts:

$$S_n=e^{n\epsilon(L+\rho/2)}+e^{n\epsilon(L-\rho/2)},\qquad
B_n=S_nU_n,\qquad P_3=B_3+CT_3=(\log Z)_3+F_3+O(\epsilon),$$
$$CT_3=-\frac{2b_0}{\epsilon}B_2+
\left(\frac{b_0^2}{\epsilon^2}-\frac{b_1}{2\epsilon}\right)B_1,
\qquad b_0=\frac{11C_A-4T}{12},\quad
b_1=\frac{34C_A^2-20C_AT-12C_FT}{48}.$$

Each lower-loop $B_n$ keeps its own $S_n$. The complete lower-loop input series are [one loop](../../benchmarks/current/results/1loop.json) and [two loops](../../benchmarks/current/results/2loop.json).

## Function definitions

Write $\zeta_n=\zeta(n)$ and suppress the argument $z$ on $M_{k,j}$ and $N_{k,j}$ in the formulas. Their raw-moment definitions are

$$r=\frac{z}{1+z},\qquad W_\beta=\log[\beta(1-\beta)],\qquad
\ell_\beta=\rho+\log\frac{1-\beta}{\beta},$$
$$M_{k,j}(z)=\int_0^1\frac{W_\beta^k\ell_\beta^j}{\beta-r}\,d\beta,
\qquad N_{k,j}(z)=M_{k,j+1}(z)-\rho M_{k,j}(z),\qquad j\ge1.$$

The numerator vanishes at $\beta=r$ and defines the removable value of the integrand. Powers are unnormalized: the definition includes no factorials. In particular,

$$M_{0,1}=-\frac{\rho^2+\pi^2}{2},\qquad M_{0,3}=-M_{0,1}^2,\qquad \pi^2=6\zeta_2.$$

The stored common-$z$ HLOG representation uses $G(a,\mathbf w;z)=\int_0^zG(\mathbf w;t)\,dt/(t-a)$, $G(\varnothing;z)=1$, with letters $0,-1$ and the standard logarithmic endpoint convention $G(0;z)=\log z$. The current four-color totals close on $M,N,\zeta_n$ and the physical scale $L$.

## Finite functions

Let $f_{3,c}(z)=F_{3,c}(0,z)$. The expressions below separate the finite shape $f_{3,c}$ from its full scale dependence.

### Finite AA: $C_FC_A^2$

$$
\begin{aligned}
F_{3,AA}(L,z)={}&f_{3,AA}(z) \\
&- \frac{121 L^{4}}{864} \\
&+ L^{3}\left(\frac{11 \zeta_{2}}{36} - \frac{445}{324}\right) \\
&+ L^{2}\left(- \frac{121 M_{0,1}}{288} - \frac{11 \zeta_{2}^{2}}{20} + \frac{899 \zeta_{2}}{288} + \frac{11 \zeta_{3}}{4} - \frac{15503}{2592}\right) \\
&+ L\left(\frac{11 M_{0,1} \zeta_{2}}{24} - \frac{445 M_{0,1}}{216} - \frac{121 M_{1,1}}{144} - \frac{11 \zeta_{2}^{2}}{5} - \frac{11 \zeta_{2} \zeta_{3}}{6} + \frac{5603 \zeta_{2}}{648} + \frac{1541 \zeta_{3}}{108} - 6 \zeta_{5} - \frac{297029}{23328}\right) .
\end{aligned}
$$

$$
\begin{aligned}
f_{3,AA}(z)={}&\frac{121 M_{0,1}^{2}}{3456} - \frac{15503 M_{0,1}}{5184} - \frac{445 M_{1,1}}{216} + \frac{11 M_{1,3}}{144} - \frac{5211949}{419904} \\
&\quad - \frac{121 M_{2,1}}{576} + \frac{11 N_{1,2}}{96} - \frac{407 \zeta_{2}^{3}}{7560} - \frac{25709 \zeta_{2}^{2}}{17280} + \frac{455435 \zeta_{2}}{46656} \\
&\quad - \frac{11 M_{0,1} \zeta_{2}^{2}}{40} + \frac{415 M_{0,1} \zeta_{2}}{576} - \frac{161 \zeta_{3}^{2}}{36} + \frac{21763 \zeta_{3}}{1944} - \frac{11 \zeta_{5}}{9} \\
&\quad + \frac{11 M_{0,1} \zeta_{3}}{4} - \frac{77 \zeta_{2} \zeta_{3}}{18} .
\end{aligned}
$$

### Finite Af: $C_FC_AT$

$$
\begin{aligned}
F_{3,Af}(L,z)={}&f_{3,Af}(z) \\
&+ \frac{11 L^{4}}{108} \\
&+ L^{3}\left(\frac{289}{324} - \frac{\zeta_{2}}{9}\right) \\
&+ L^{2}\left(\frac{11 M_{0,1}}{36} - \frac{53 \zeta_{2}}{36} + \frac{2051}{648}\right) \\
&+ L\left(- \frac{M_{0,1} \zeta_{2}}{6} + \frac{289 M_{0,1}}{216} + \frac{11 M_{1,1}}{18} + \frac{6 \zeta_{2}^{2}}{5} - \frac{3013 \zeta_{2}}{648} - \frac{113 \zeta_{3}}{54} + \frac{31313}{5832}\right) .
\end{aligned}
$$

$$
\begin{aligned}
f_{3,Af}(z)={}&- \frac{11 M_{0,1}^{2}}{432} + \frac{2051 M_{0,1}}{1296} + \frac{289 M_{1,1}}{216} - \frac{M_{1,3}}{36} + \frac{412765}{104976} \\
&\quad + \frac{11 M_{2,1}}{72} - \frac{N_{1,2}}{24} + \frac{2239 \zeta_{2}^{2}}{2160} - \frac{56071 \zeta_{2}}{11664} - \frac{76 \zeta_{3}}{81} \\
&\quad - \frac{M_{0,1} \zeta_{2}}{8} - \frac{M_{0,1} \zeta_{3}}{2} + \frac{7 \zeta_{2} \zeta_{3}}{18} + 2 \zeta_{5} .
\end{aligned}
$$

### Finite Ff: $C_F^2T$

$$
\begin{aligned}
F_{3,Ff}(L,z)={}&f_{3,Ff}(z) \\
&+ \frac{L^{3}}{12} \\
&+ L^{2}\left(\frac{55}{48} - \zeta_{3}\right) \\
&+ L\left(\frac{M_{0,1}}{8} - \frac{2 \zeta_{2}^{2}}{5} - \frac{3 \zeta_{2}}{8} - \frac{19 \zeta_{3}}{9} + \frac{1711}{432}\right) .
\end{aligned}
$$

$$
\begin{aligned}
f_{3,Ff}(z)={}&\frac{55 M_{0,1}}{96} + \frac{M_{1,1}}{8} - \frac{19 \zeta_{2}^{2}}{45} - \frac{385 \zeta_{2}}{288} + \frac{42727}{7776} \\
&\quad - \frac{M_{0,1} \zeta_{3}}{2} + \frac{7 \zeta_{2} \zeta_{3}}{6} - \frac{355 \zeta_{3}}{162} - \frac{14 \zeta_{5}}{9} .
\end{aligned}
$$

### Finite ff: $C_FT^2$

$$
\begin{aligned}
F_{3,ff}(L,z)={}&f_{3,ff}(z) \\
&- \frac{L^{4}}{54} \\
&- \frac{10 L^{3}}{81} \\
&+ L^{2}\left(- \frac{M_{0,1}}{18} + \frac{\zeta_{2}}{6} - \frac{25}{81}\right) \\
&+ L\left(- \frac{5 M_{0,1}}{27} - \frac{M_{1,1}}{9} + \frac{5 \zeta_{2}}{9} - \frac{4 \zeta_{3}}{9} - \frac{232}{729}\right) .
\end{aligned}
$$

$$
\begin{aligned}
f_{3,ff}(z)={}&\frac{M_{0,1}^{2}}{216} - \frac{25 M_{0,1}}{162} - \frac{5 M_{1,1}}{27} - \frac{M_{2,1}}{36} + \frac{32}{6561} \\
&\quad - \frac{M_{0,1} \zeta_{2}}{36} - \frac{149 \zeta_{2}^{2}}{1080} + \frac{77 \zeta_{2}}{162} - \frac{110 \zeta_{3}}{243} .
\end{aligned}
$$

## All-Laurent coefficients

For each color channel $c\in\{AA,Af,Ff,ff\}$ define

$$U_{3,c}=\sum_{p=-6}^{0}u_{c,p}\epsilon^p+O(\epsilon),\quad
B_{3,c}=\sum_{p=-6}^{0}b_{c,p}\epsilon^p+O(\epsilon),\quad
P_{3,c}=\sum_{p=-6}^{-1}q_{c,p}\epsilon^p+F_{3,c}+O(\epsilon).$$

Every coefficient, including zero entries, is displayed below. The counterterm coefficient is $[\epsilon^p]CT_{3,c}=q_{c,p}-b_{c,p}$ for $p<0$ and $F_{3,c}-b_{c,0}$ for $p=0$; its complete evaluated values are also stored in the exact data.

### Laurent AA: $C_FC_A^2$

**Scale-stripped sum.**

$$
\begin{aligned}
u_{AA,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{AA,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{AA,-4}={}&- \frac{121}{1296} .
\end{aligned}
$$

$$
\begin{aligned}
u_{AA,-3}={}&\frac{11 \zeta_{2}}{108} - \frac{1627}{3888} .
\end{aligned}
$$

$$
\begin{aligned}
u_{AA,-2}={}&- \frac{121 M_{0,1}}{288} - \frac{11 \zeta_{2}^{2}}{180} + \frac{173 \zeta_{2}}{2592} + \frac{55 \zeta_{3}}{54} - \frac{11093}{7776} .
\end{aligned}
$$

$$
\begin{aligned}
u_{AA,-1}={}&- \frac{1627 M_{0,1}}{864} - \frac{121 M_{1,1}}{192} + \frac{77 \zeta_{2}^{2}}{240} - \frac{1685 \zeta_{2}}{7776} - \frac{617525}{139968} \\
&\quad + \frac{11 M_{0,1} \zeta_{2}}{24} - \frac{11 \zeta_{2} \zeta_{3}}{36} + \frac{329 \zeta_{3}}{144} - \zeta_{5} .
\end{aligned}
$$

$$
\begin{aligned}
u_{AA,0}={}&- \frac{605 M_{0,1}^{2}}{768} - \frac{11093 M_{0,1}}{1728} - \frac{1627 M_{1,1}}{576} + \frac{11 M_{1,3}}{96} - \frac{10984045}{839808} \\
&\quad - \frac{121 M_{2,1}}{256} + \frac{11 N_{1,2}}{64} - \frac{1661 \zeta_{2}^{3}}{3780} + \frac{12173 \zeta_{2}^{2}}{17280} - \frac{61565 \zeta_{2}}{46656} \\
&\quad - \frac{11 M_{0,1} \zeta_{2}^{2}}{40} - \frac{2921 M_{0,1} \zeta_{2}}{1152} - \frac{161 \zeta_{3}^{2}}{72} + \frac{19093 \zeta_{3}}{3888} + \frac{671 \zeta_{5}}{144} \\
&\quad + \frac{319 M_{0,1} \zeta_{3}}{48} + \frac{11 \zeta_{2} \zeta_{3}}{9} .
\end{aligned}
$$

**Scale-restored bare exponent.**

$$
\begin{aligned}
b_{AA,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{AA,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{AA,-4}={}&- \frac{121}{648} .
\end{aligned}
$$

$$
\begin{aligned}
b_{AA,-3}={}&- \frac{121 L}{216} + \frac{11 \zeta_{2}}{54} - \frac{1627}{1944} .
\end{aligned}
$$

$$
\begin{aligned}
b_{AA,-2}={}&- \frac{121 L^{2}}{144} - \frac{1627 L}{648} - \frac{121 M_{0,1}}{288} + \frac{3613 \zeta_{2}}{2592} - \frac{11093}{3888} \\
&\quad + \frac{11 L \zeta_{2}}{18} - \frac{11 \zeta_{2}^{2}}{90} + \frac{55 \zeta_{3}}{27} .
\end{aligned}
$$

$$
\begin{aligned}
b_{AA,-1}={}&- \frac{121 L^{3}}{144} - \frac{1627 L^{2}}{432} - \frac{11093 L}{1296} - \frac{1627 M_{0,1}}{864} - \frac{617525}{69984} \\
&\quad - \frac{121 M_{1,1}}{96} - \frac{11 \zeta_{2}^{2}}{15} + \frac{40559 \zeta_{2}}{7776} + \frac{329 \zeta_{3}}{72} - 2 \zeta_{5} \\
&\quad + \frac{11 L^{2} \zeta_{2}}{12} - \frac{121 L M_{0,1}}{96} - \frac{11 L \zeta_{2}^{2}}{30} + \frac{3613 L \zeta_{2}}{864} + \frac{55 L \zeta_{3}}{9} \\
&\quad + \frac{11 M_{0,1} \zeta_{2}}{24} - \frac{11 \zeta_{2} \zeta_{3}}{18} .
\end{aligned}
$$

$$
\begin{aligned}
b_{AA,0}={}&- \frac{121 L^{4}}{192} - \frac{1627 L^{3}}{432} - \frac{11093 L^{2}}{864} - \frac{617525 L}{23328} - \frac{10984045}{419904} \\
&\quad + \frac{121 M_{0,1}^{2}}{768} - \frac{11093 M_{0,1}}{1728} - \frac{1627 M_{1,1}}{288} + \frac{11 M_{1,3}}{48} - \frac{121 M_{2,1}}{128} \\
&\quad + \frac{11 N_{1,2}}{32} - \frac{407 \zeta_{2}^{3}}{7560} - \frac{31453 \zeta_{2}^{2}}{34560} + \frac{775403 \zeta_{2}}{46656} + \frac{19093 \zeta_{3}}{1944} \\
&\quad - \frac{1627 L M_{0,1}}{288} - \frac{121 L M_{1,1}}{32} + \frac{40559 L \zeta_{2}}{2592} - \frac{161 \zeta_{3}^{2}}{36} + \frac{671 \zeta_{5}}{72} \\
&\quad - \frac{121 L^{2} M_{0,1}}{64} + \frac{3613 L^{2} \zeta_{2}}{576} - \frac{11 L \zeta_{2}^{2}}{5} + \frac{329 L \zeta_{3}}{24} - 6 L \zeta_{5} \\
&\quad + \frac{11 L^{3} \zeta_{2}}{12} - \frac{11 L^{2} \zeta_{2}^{2}}{20} + \frac{55 L^{2} \zeta_{3}}{6} - \frac{11 M_{0,1} \zeta_{2}^{2}}{40} - \frac{743 M_{0,1} \zeta_{2}}{1152} \\
&\quad + \frac{11 L M_{0,1} \zeta_{2}}{8} - \frac{11 L \zeta_{2} \zeta_{3}}{6} + \frac{209 M_{0,1} \zeta_{3}}{24} - \frac{407 \zeta_{2} \zeta_{3}}{36} .
\end{aligned}
$$

**Pole coefficients after coupling renormalization.**

$$
\begin{aligned}
q_{AA,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{AA,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{AA,-4}={}&- \frac{1331}{1296} .
\end{aligned}
$$

$$
\begin{aligned}
q_{AA,-3}={}&- \frac{121 L}{216} - \frac{55 \zeta_{2}}{216} + \frac{6133}{3888} .
\end{aligned}
$$

$$
\begin{aligned}
q_{AA,-2}={}&\frac{1043 L}{648} - \frac{11 \zeta_{2}^{2}}{90} + \frac{173 \zeta_{2}}{1296} - \frac{253 \zeta_{3}}{216} + \frac{2239}{3888} \\
&\quad - \frac{11 L \zeta_{2}}{36} .
\end{aligned}
$$

$$
\begin{aligned}
q_{AA,-1}={}&- \frac{245 L}{144} - \frac{11 \zeta_{2}^{2}}{15} + \frac{6325 \zeta_{2}}{3888} + \frac{329 \zeta_{3}}{72} - \frac{136781}{69984} \\
&\quad - \frac{11 L \zeta_{2}^{2}}{30} + \frac{67 L \zeta_{2}}{54} - \frac{11 L \zeta_{3}}{36} - \frac{11 \zeta_{2} \zeta_{3}}{18} - 2 \zeta_{5} .
\end{aligned}
$$

The zero-order coefficient is $q_{AA,0}=F_{3,AA}$, given above.

### Laurent Af: $C_FC_AT$

**Scale-stripped sum.**

$$
\begin{aligned}
u_{Af,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{Af,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{Af,-4}={}&\frac{11}{162} .
\end{aligned}
$$

$$
\begin{aligned}
u_{Af,-3}={}&\frac{533}{1944} - \frac{\zeta_{2}}{27} .
\end{aligned}
$$

$$
\begin{aligned}
u_{Af,-2}={}&\frac{11 M_{0,1}}{36} + \frac{13 \zeta_{2}}{324} - \frac{7 \zeta_{3}}{27} + \frac{1633}{1944} .
\end{aligned}
$$

$$
\begin{aligned}
u_{Af,-1}={}&\frac{533 M_{0,1}}{432} + \frac{11 M_{1,1}}{24} - \frac{\zeta_{2}^{2}}{20} + \frac{1187 \zeta_{2}}{3888} + \frac{82097}{34992} \\
&\quad - \frac{M_{0,1} \zeta_{2}}{6} - \frac{91 \zeta_{3}}{324} .
\end{aligned}
$$

$$
\begin{aligned}
u_{Af,0}={}&\frac{55 M_{0,1}^{2}}{96} + \frac{1633 M_{0,1}}{432} + \frac{533 M_{1,1}}{288} - \frac{M_{1,3}}{24} + \frac{1323997}{209952} \\
&\quad + \frac{11 M_{2,1}}{32} - \frac{N_{1,2}}{16} + \frac{137 \zeta_{2}^{2}}{2160} + \frac{12469 \zeta_{2}}{11664} - \frac{5 \zeta_{3}}{216} \\
&\quad + \frac{323 M_{0,1} \zeta_{2}}{144} - \frac{23 M_{0,1} \zeta_{3}}{12} - \frac{5 \zeta_{2} \zeta_{3}}{18} - \frac{11 \zeta_{5}}{12} .
\end{aligned}
$$

**Scale-restored bare exponent.**

$$
\begin{aligned}
b_{Af,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{Af,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{Af,-4}={}&\frac{11}{81} .
\end{aligned}
$$

$$
\begin{aligned}
b_{Af,-3}={}&\frac{11 L}{27} - \frac{2 \zeta_{2}}{27} + \frac{533}{972} .
\end{aligned}
$$

$$
\begin{aligned}
b_{Af,-2}={}&\frac{11 L^{2}}{18} + \frac{533 L}{324} + \frac{11 M_{0,1}}{36} - \frac{271 \zeta_{2}}{324} + \frac{1633}{972} \\
&\quad - \frac{2 L \zeta_{2}}{9} - \frac{14 \zeta_{3}}{27} .
\end{aligned}
$$

$$
\begin{aligned}
b_{Af,-1}={}&\frac{11 L^{3}}{18} + \frac{533 L^{2}}{216} + \frac{1633 L}{324} + \frac{533 M_{0,1}}{432} + \frac{82097}{17496} \\
&\quad + \frac{11 L M_{0,1}}{12} + \frac{11 M_{1,1}}{12} + \frac{2 \zeta_{2}^{2}}{5} - \frac{12017 \zeta_{2}}{3888} - \frac{91 \zeta_{3}}{162} \\
&\quad - \frac{L^{2} \zeta_{2}}{3} - \frac{271 L \zeta_{2}}{108} - \frac{14 L \zeta_{3}}{9} - \frac{M_{0,1} \zeta_{2}}{6} .
\end{aligned}
$$

$$
\begin{aligned}
b_{Af,0}={}&\frac{11 L^{4}}{24} + \frac{533 L^{3}}{216} + \frac{1633 L^{2}}{216} + \frac{82097 L}{5832} + \frac{1323997}{104976} \\
&\quad - \frac{11 M_{0,1}^{2}}{96} + \frac{1633 M_{0,1}}{432} + \frac{533 M_{1,1}}{144} - \frac{M_{1,3}}{12} + \frac{11 M_{2,1}}{16} \\
&\quad - \frac{N_{1,2}}{8} + \frac{2663 \zeta_{2}^{2}}{4320} - \frac{107335 \zeta_{2}}{11664} - \frac{5 \zeta_{3}}{108} - \frac{11 \zeta_{5}}{6} \\
&\quad + \frac{533 L M_{0,1}}{144} + \frac{11 L M_{1,1}}{4} + \frac{6 L \zeta_{2}^{2}}{5} - \frac{12017 L \zeta_{2}}{1296} - \frac{91 L \zeta_{3}}{54} \\
&\quad - \frac{L^{3} \zeta_{2}}{3} + \frac{11 L^{2} M_{0,1}}{8} - \frac{271 L^{2} \zeta_{2}}{72} - \frac{7 L^{2} \zeta_{3}}{3} + \frac{125 M_{0,1} \zeta_{2}}{144} \\
&\quad - \frac{L M_{0,1} \zeta_{2}}{2} - \frac{8 M_{0,1} \zeta_{3}}{3} + \frac{53 \zeta_{2} \zeta_{3}}{18} .
\end{aligned}
$$

**Pole coefficients after coupling renormalization.**

$$
\begin{aligned}
q_{Af,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{Af,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{Af,-4}={}&\frac{121}{162} .
\end{aligned}
$$

$$
\begin{aligned}
q_{Af,-3}={}&\frac{11 L}{27} + \frac{5 \zeta_{2}}{54} - \frac{485}{486} .
\end{aligned}
$$

$$
\begin{aligned}
q_{Af,-2}={}&\frac{L \zeta_{2}}{9} - \frac{167 L}{162} + \frac{13 \zeta_{2}}{162} + \frac{35 \zeta_{3}}{54} - \frac{503}{972} .
\end{aligned}
$$

$$
\begin{aligned}
q_{Af,-1}={}&\frac{209 L}{324} + \frac{2 \zeta_{2}^{2}}{5} - \frac{707 \zeta_{2}}{972} - \frac{91 \zeta_{3}}{162} + \frac{5921}{17496} \\
&\quad - \frac{10 L \zeta_{2}}{27} + \frac{7 L \zeta_{3}}{9} .
\end{aligned}
$$

The zero-order coefficient is $q_{Af,0}=F_{3,Af}$, given above.

### Laurent Ff: $C_F^2T$

**Scale-stripped sum.**

$$
\begin{aligned}
u_{Ff,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{Ff,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{Ff,-4}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{Ff,-3}={}&\frac{1}{72} .
\end{aligned}
$$

$$
\begin{aligned}
u_{Ff,-2}={}&\frac{55}{432} - \frac{\zeta_{3}}{9} .
\end{aligned}
$$

$$
\begin{aligned}
u_{Ff,-1}={}&\frac{M_{0,1}}{16} - \frac{\zeta_{2}^{2}}{15} + \frac{\zeta_{2}}{48} - \frac{19 \zeta_{3}}{54} + \frac{1711}{2592} .
\end{aligned}
$$

$$
\begin{aligned}
u_{Ff,0}={}&\frac{55 M_{0,1}}{96} + \frac{3 M_{1,1}}{32} - \frac{19 \zeta_{2}^{2}}{90} + \frac{55 \zeta_{2}}{288} + \frac{42727}{15552} \\
&\quad - \frac{M_{0,1} \zeta_{3}}{2} - \frac{\zeta_{2} \zeta_{3}}{6} - \frac{683 \zeta_{3}}{648} - \frac{7 \zeta_{5}}{9} .
\end{aligned}
$$

**Scale-restored bare exponent.**

$$
\begin{aligned}
b_{Ff,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{Ff,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{Ff,-4}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{Ff,-3}={}&\frac{1}{36} .
\end{aligned}
$$

$$
\begin{aligned}
b_{Ff,-2}={}&\frac{L}{12} - \frac{2 \zeta_{3}}{9} + \frac{55}{216} .
\end{aligned}
$$

$$
\begin{aligned}
b_{Ff,-1}={}&\frac{L^{2}}{8} + \frac{55 L}{72} + \frac{M_{0,1}}{16} - \frac{7 \zeta_{2}}{48} + \frac{1711}{1296} \\
&\quad - \frac{2 L \zeta_{3}}{3} - \frac{2 \zeta_{2}^{2}}{15} - \frac{19 \zeta_{3}}{27} .
\end{aligned}
$$

$$
\begin{aligned}
b_{Ff,0}={}&\frac{L^{3}}{8} + \frac{55 L^{2}}{48} + \frac{1711 L}{432} + \frac{55 M_{0,1}}{96} + \frac{42727}{7776} \\
&\quad + \frac{3 M_{1,1}}{16} - \frac{19 \zeta_{2}^{2}}{45} - \frac{385 \zeta_{2}}{288} - \frac{683 \zeta_{3}}{324} - \frac{14 \zeta_{5}}{9} \\
&\quad - L^{2} \zeta_{3} + \frac{3 L M_{0,1}}{16} - \frac{2 L \zeta_{2}^{2}}{5} - \frac{7 L \zeta_{2}}{16} - \frac{19 L \zeta_{3}}{9} \\
&\quad - \frac{M_{0,1} \zeta_{3}}{2} + \frac{7 \zeta_{2} \zeta_{3}}{6} .
\end{aligned}
$$

**Pole coefficients after coupling renormalization.**

$$
\begin{aligned}
q_{Ff,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{Ff,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{Ff,-4}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{Ff,-3}={}&- \frac{2}{9} .
\end{aligned}
$$

$$
\begin{aligned}
q_{Ff,-2}={}&- \frac{L}{6} - \frac{2 \zeta_{3}}{9} + \frac{55}{216} .
\end{aligned}
$$

$$
\begin{aligned}
q_{Ff,-1}={}&\frac{55 L}{72} - \frac{2 \zeta_{2}^{2}}{15} - \frac{\zeta_{2}}{12} - \frac{19 \zeta_{3}}{27} + \frac{1711}{1296} \\
&\quad - \frac{2 L \zeta_{3}}{3} .
\end{aligned}
$$

The zero-order coefficient is $q_{Ff,0}=F_{3,Ff}$, given above.

### Laurent ff: $C_FT^2$

**Scale-stripped sum.**

$$
\begin{aligned}
u_{ff,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{ff,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
u_{ff,-4}={}&- \frac{1}{81} .
\end{aligned}
$$

$$
\begin{aligned}
u_{ff,-3}={}&- \frac{10}{243} .
\end{aligned}
$$

$$
\begin{aligned}
u_{ff,-2}={}&- \frac{M_{0,1}}{18} - \frac{\zeta_{2}}{54} - \frac{1}{9} .
\end{aligned}
$$

$$
\begin{aligned}
u_{ff,-1}={}&- \frac{5 M_{0,1}}{27} - \frac{M_{1,1}}{12} - \frac{5 \zeta_{2}}{81} - \frac{7 \zeta_{3}}{81} - \frac{608}{2187} .
\end{aligned}
$$

$$
\begin{aligned}
u_{ff,0}={}&- \frac{5 M_{0,1}^{2}}{48} - \frac{M_{0,1}}{2} - \frac{5 M_{1,1}}{18} - \frac{M_{2,1}}{16} - \frac{4376}{6561} \\
&\quad - \frac{11 M_{0,1} \zeta_{2}}{24} - \frac{67 \zeta_{2}^{2}}{1080} - \frac{\zeta_{2}}{6} - \frac{70 \zeta_{3}}{243} .
\end{aligned}
$$

**Scale-restored bare exponent.**

$$
\begin{aligned}
b_{ff,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{ff,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
b_{ff,-4}={}&- \frac{2}{81} .
\end{aligned}
$$

$$
\begin{aligned}
b_{ff,-3}={}&- \frac{2 L}{27} - \frac{20}{243} .
\end{aligned}
$$

$$
\begin{aligned}
b_{ff,-2}={}&- \frac{L^{2}}{9} - \frac{20 L}{81} - \frac{M_{0,1}}{18} + \frac{7 \zeta_{2}}{54} - \frac{2}{9} .
\end{aligned}
$$

$$
\begin{aligned}
b_{ff,-1}={}&- \frac{L^{3}}{9} - \frac{10 L^{2}}{27} - \frac{2 L}{3} - \frac{5 M_{0,1}}{27} - \frac{1216}{2187} \\
&\quad - \frac{L M_{0,1}}{6} + \frac{7 L \zeta_{2}}{18} - \frac{M_{1,1}}{6} + \frac{35 \zeta_{2}}{81} - \frac{14 \zeta_{3}}{81} .
\end{aligned}
$$

$$
\begin{aligned}
b_{ff,0}={}&- \frac{L^{4}}{12} - \frac{10 L^{3}}{27} - L^{2} - \frac{1216 L}{729} - \frac{8752}{6561} \\
&\quad + \frac{M_{0,1}^{2}}{48} - \frac{M_{0,1}}{2} - \frac{5 M_{1,1}}{9} - \frac{M_{2,1}}{8} + \frac{7 \zeta_{2}}{6} \\
&\quad - \frac{5 L M_{0,1}}{9} - \frac{L M_{1,1}}{2} + \frac{35 L \zeta_{2}}{27} - \frac{133 \zeta_{2}^{2}}{2160} - \frac{140 \zeta_{3}}{243} \\
&\quad - \frac{L^{2} M_{0,1}}{4} + \frac{7 L^{2} \zeta_{2}}{12} - \frac{14 L \zeta_{3}}{27} - \frac{5 M_{0,1} \zeta_{2}}{24} .
\end{aligned}
$$

**Pole coefficients after coupling renormalization.**

$$
\begin{aligned}
q_{ff,-6}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{ff,-5}={}&0 .
\end{aligned}
$$

$$
\begin{aligned}
q_{ff,-4}={}&- \frac{11}{81} .
\end{aligned}
$$

$$
\begin{aligned}
q_{ff,-3}={}&\frac{25}{243} - \frac{2 L}{27} .
\end{aligned}
$$

$$
\begin{aligned}
q_{ff,-2}={}&\frac{10 L}{81} - \frac{\zeta_{2}}{27} + \frac{10}{81} .
\end{aligned}
$$

$$
\begin{aligned}
q_{ff,-1}={}&\frac{2 L}{81} + \frac{5 \zeta_{2}}{81} - \frac{14 \zeta_{3}}{81} + \frac{260}{2187} .
\end{aligned}
$$

The zero-order coefficient is $q_{ff,0}=F_{3,ff}$, given above.

## Renormalization-group coefficients

The cusp and noncusp coefficients have the same four-color decomposition, $k_3=C_F\sum_c\mathcal C_ck_{3,c}$ and $h_3=C_F\sum_c\mathcal C_ch_{3,c}$. In the convention above,

$$\partial_LF_3-2b_0F_2-b_1F_1=-2k_3L-h_3.$$

**AA.**

$$k_{3,AA}=\frac{11 \zeta_{2}^{2}}{20} - \frac{67 \zeta_{2}}{36} + \frac{11 \zeta_{3}}{24} + \frac{245}{96},\qquad h_{3,AA}=\frac{11 \zeta_{2}^{2}}{5} + \frac{11 \zeta_{2} \zeta_{3}}{6} - \frac{6325 \zeta_{2}}{1296} - \frac{329 \zeta_{3}}{24} + 6 \zeta_{5} + \frac{136781}{23328}.$$

**Af.**

$$k_{3,Af}=\frac{5 \zeta_{2}}{9} - \frac{7 \zeta_{3}}{6} - \frac{209}{216},\qquad h_{3,Af}=- \frac{6 \zeta_{2}^{2}}{5} + \frac{707 \zeta_{2}}{324} + \frac{91 \zeta_{3}}{54} - \frac{5921}{5832}.$$

**Ff.**

$$k_{3,Ff}=\zeta_{3} - \frac{55}{48},\qquad h_{3,Ff}=\frac{2 \zeta_{2}^{2}}{5} + \frac{\zeta_{2}}{4} + \frac{19 \zeta_{3}}{9} - \frac{1711}{432}.$$

**ff.**

$$k_{3,ff}=- \frac{1}{27},\qquad h_{3,ff}=- \frac{5 \zeta_{2}}{27} + \frac{14 \zeta_{3}}{27} - \frac{260}{729}.$$

## Exact data and calculation record

The [three-loop JSON](../../benchmarks/current/results/3loop.json) stores these formulas under `channels.AA`, `channels.Af`, `channels.Ff` and `channels.ff`: `U`, `bare`, `counterterm`, `post_counterterm`, `logZ`, `F`, `f_at_L0`, `k3` and `h3`. `native_common_z_HLOG_channels` supplies the corresponding HLOG construction. The [current source catalogue](../../benchmarks/current/catalog.json) connects all163 source orbits to their inputs and selected integral objects.

The [source/joint ledger](../2026-11-02-current/source-ledger.json) and [matter ledger](../2026-11-02-current/matter-ledger.json) record their consumption. The accepted integral derivations, exact function dictionary, source assembly and RG checks form the calculation record. Read the [integral derivations](../../integrals/README.md), [calculation validation](../../validation/README.md) and [formula-presentation check](../../validation/results-presentation-20261103/README.md) for the corresponding steps.

This article is generated from the current exact data by [the formula renderer](../../tools/render_four_color_results.py). Its check reproduces every displayed coefficient, the finite slices, scale restoration, lower-loop counterterms, pole separation and the four-color finite RG relation.

```text
python -B tools/render_four_color_results.py --check
python -m workflow.current_calculation result 3
```

Earlier releases and dated research records are indexed in the [historical archive](../../docs/archive/README.md).
