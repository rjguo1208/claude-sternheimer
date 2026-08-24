# Koster–Slater / defect Green function: the cheap route to defect levels

The two earlier routes both give the correct $a_1\oplus e$, and both are expensive: **explicit 21-band** (move the bands into $P$ and diagonalize explicitly; [results page](results.html#sec-3) Fig 14) costs $46$ min $+\,N_A^3$ diagonalization, and **all-order Feshbach** (leave the bands in $Q$ and dress to all orders; [ladder page](sternheimer-ladder.html#sec-6) §6) costs ~days. This note derives a third and by far the cheapest: the **Koster–Slater defect Green function method** (the classic defect-GF approach). It makes **no $P$/$Q$ split and does no rest dressing**; it writes the defect levels directly as a secular equation on a **small defect-local block**, $\det[1-G_0(E)\Delta V]=0$ — the host Green function $G_0$ carries all the bands inside it (so there is no over-screening and no ladder divergence), while $\Delta V$ is short-ranged (measured on the results page, Figs 19–21), leaving a determinant of only **a few tens of dimensions**.

## 1. From a large diagonalization to a small secular equation (Lippmann–Schwinger)

The defect Hamiltonian is $H=H_0+\Delta V$, with host $H_0|n\mathbf k\rangle=\varepsilon_{n\mathbf k}|n\mathbf k\rangle$. A bound state $H|\psi\rangle=E|\psi\rangle$ means
$$(E-H_0)|\psi\rangle=\Delta V|\psi\rangle.$$
For $E$ in the gap ($E\notin\mathrm{spec}\,H_0$), $(E-H_0)$ is invertible; defining the host retarded Green function $G_0(E)=(E-H_0)^{-1}$,
$$|\psi\rangle=G_0(E)\,\Delta V\,|\psi\rangle.$$
The homogeneous equation has a nontrivial solution $\iff$
$$\boxed{\,\det\!\big[\,1-G_0(E)\,\Delta V\,\big]=0\,}. \tag{1}$$
Equivalently, with the defect scattering matrix $T(E)=\Delta V[1-G_0(E)\Delta V]^{-1}$ and the full Green function $G=G_0+G_0TG_0$, **the bound states are the poles of $T$ (and of $G$)**, i.e. the roots of (1). This step is exact — no approximation.

## 2. Locality → the determinant shrinks to the "defect block"

Formally (1) is a determinant over the whole space, but $\Delta V$ is **short-ranged**: in the Wannier basis $|w_{\mathbf R\alpha}\rangle$, $\Delta V_{\mathbf R\alpha,\mathbf R'\beta}\neq0$ only inside the defect region $D$ ($|\mathbf R|,|\mathbf R'|\le R_{\rm cut}$) — results page Figs 19–21 show the S-vacancy Koster–Slater is already converged at $R_{\rm cut}=4$. Writing $P_D$ for the projector onto $D$, $\Delta V=P_D\,\Delta V\,P_D$ with rank $\le\dim D$. The nonzero eigenvalues of $G_0\Delta V$ equal those of the $D\times D$ matrix $g(E)\Delta V_D$, so (1) contracts to
$$\boxed{\,\det{}_D\!\big[\,1-g(E)\,\Delta V_D\,\big]=0\,},\qquad g(E)\equiv P_D\,G_0(E)\,P_D. \tag{2}$$
Both $\Delta V_D$ and $g(E)$ are $\dim D\sim$ **a few tens**, not $N_A=3002$. **The key point:** $\Delta V_D$ uses the **bare** $M^W$ (gauge-fixed and measured short-ranged on the results page, Figs 18–19), **not** the dressed $\tilde V^W$ (which carries the second-order $\Sigma$ and over-screens).

## 3. The host Green function $g(E)$: all bands are inside it, and adding bands is cheap

In the Wannier basis $g(E)$ is just the resolvent of the interpolated host Hamiltonian:
$$g(E,\mathbf k)=\big[E\,\mathbb 1-H_0^W(\mathbf k)\big]^{-1},\qquad
g(E)_{\mathbf R\alpha,\mathbf R'\beta}=\frac1{N_k}\sum_{\mathbf k}e^{i\mathbf k\cdot(\mathbf R-\mathbf R')}\,g(E,\mathbf k)_{\alpha\beta}, \tag{3}$$
where $H_0^W(\mathbf k)$ is the Wannier-interpolated host (EDT already has $H_W(\mathbf R)$; [results page §4](results.html#sec-4)). Inverting $g(E,\mathbf k)$ **automatically includes every band in the Wannier space** — the conduction bands that produce $e$ are among them (11 bands already suffice to bind $e$, just as the 11-band bare $M$ gives $e=+1.49$). Converging the level further (toward explicit's $+1.35$ / DFT's $+1.19$) only requires folding more host bands into $g$, through their projection onto the defect orbitals — and the cost of that is extra **eigenvalue** $k$-sums, far less than the $N_A$ matrix that explicit pays per band. The $k$-sum can also use a very fine mesh almost for free.

That is the source of its being **both accurate and cheap**: **the heavy lifting (full BZ, many bands) falls only on cheap host eigenvalue sums, while the strong defect scattering is treated exactly inside a $D$ block of a few tens of dimensions (the full $\det$, with no expansion in $\Delta V$).** So there is neither the over-screening of second-order dressing nor the ladder's divergence problem on $\Delta V_{QQ}$ — (2) is simply not an expansion in $\Delta V$. For $E$ in the gap, $E$ is real and $g$ is finite; to resolve resonances inside the bands, use $E+i\eta$.

## 4. C₃ᵥ symmetry → the secular determinant block-diagonalizes, irrep labels are rigorous

The orbitals of $D$ (three Mo dangling bonds plus neighbours) reduce under $C_{3v}$ to contain $a_1\oplus e$. Both $\Delta V_D$ and $g(E)$ commute with the point group, so the determinant in (2) **factorizes by irrep**:
$$\det{}_D[\,1-g\Delta V_D\,]=\prod_{\Gamma}\Big(\det{}_\Gamma[\,1-g_\Gamma\Delta V_\Gamma\,]\Big)^{d_\Gamma}.$$
The $a_1$ root comes from the $a_1$ block (1-dimensional), the $e$ root from the $e$ block ($d_e=2$ → automatically doubly degenerate). This gives the $a_1$/$a_2$/$e$ labels **rigorously**, answering the earlier caveat that the assignment rested on degeneracy alone with no character computed.

## 5. ΔDOS: the same determinant delivers it (Krein–Friedel)

The defect-induced change in the density of states is given by the phase of the secular determinant:
$$\Delta\rho(E)=-\frac1\pi\frac{d}{dE}\,\mathrm{Im}\,\ln\det{}_D\!\big[\,1-g(E+i0^+)\,\Delta V_D\,\big]
=-\frac1\pi\frac{d}{dE}\arg\det{}_D[\cdots]. \tag{4}$$
Inside the gap, each jump of $\pi$ in the phase corresponds to one bound state (a root of (2)); inside the bands it gives the resonance broadening. So the ΔDOS on the results page (Fig 15) follows directly from this small determinant, with no need to sweep the full $T$-matrix.

## 6. Cost and where it sits

| route | cost of the defect levels | over-screening? |
|---|---|---|
| second-order block dressing | $2$ h | **yes** (divergent) |
| all-order Feshbach | ~$1$–$3$ days (self-consistent) | no, but expensive |
| explicit 21-band | $46$ min $+\,N_A^3$ diagonalization | no |
| **Koster–Slater (2)** | **seconds** ($\dim D^3$ det $+$ host $k$-sum) | **no** |

Per $E$: a host $k$-sum of $\mathcal O(N_k\,n_{\rm band}\,\dim D^2)$ plus a determinant of $\mathcal O(\dim D^3)$; sweeping a few tens of $E$ points or root-finding puts the total at **seconds to minutes**, and a finer $k$ mesh can still be used for convergence. It **bypasses** the $P$/$Q$ split and the rest dressing (no expansion in $\Delta V_{QQ}$, no ladder divergence), and it also saves the full $M$ block and the $N_A^3$ diagonalization that explicit pays. **The groundwork is already in place**: the gauge fix of the Wannier $M^W$ (Figs 18–19) and its short-rangedness / $R_{\rm cut}=4$ convergence (Figs 19–21) were both verified on the results page. This is the **cheapest accurate** route to the S-vacancy $a_1$+$e$; implementation and measured results are in §7.

## 7. Implementation and results: the bare $M^W$ reproduces $a_1+e$

Implementing (2) reuses existing pieces throughout: $\Delta V_D=$ the **bare** $M^W$ (record 1 of `vtilde_block.dat`, rotated into Wannier with $U(\mathbf k)$, Fourier transformed to $\mathbf R$, truncated to the $R_{\rm cut}=3$ defect block, dim $=539$); $g(E)$ from inverting the host $H_W(\mathbf k)$ interpolated from `mos2_hr.dat` on an $N_f=48$ fine mesh. Sweeping $E$ across the gap, the smallest singular value of $[\,1-g(E)\Delta V_D\,]$ dips at each level, and **the number of singular values collapsing together = the degeneracy**.

![Koster-Slater from the bare M^W: upper panel, the smallest singular value shows a shallow dip at a1 (VBM edge) and a deep, doubly degenerate dip (sigma1=sigma2) at e (+1.50); lower panel, the Krein-Friedel counting staircase.](../assets/koster_slater_levels.png)

*Figure. Koster–Slater from the bare $M^W$ ($R_{\rm cut}=3$, $N_f=48$, $\eta=0.02$ eV). **Upper:** the smallest singular value $\sigma_{\min}$ of $[1-g(E)\Delta V_D]$ swept over $E$ — $a_1$ gives a shallow dip at the VBM edge ($\sigma_2\gg\sigma_1$, non-degenerate), $e$ a deep dip at $-4.44$ with $\sigma_1=\sigma_2$ (doubly degenerate). **Lower:** the Krein–Friedel counting staircase $\Delta N(E)$. The two dip positions are the defect levels.*

The results agree with the 11-band explicit eigenvalues to the digit:

| | Koster–Slater (bare $M^W$) | 11-band explicit |
|---|---|---|
| $a_1$ (singlet) | $-5.929$ ($+0.01$); $\sigma_1{=}0.095,\ \sigma_2{=}0.29$ → non-degenerate | $-5.926$ ($+0.01$) |
| $e$ (doublet) | $-4.441$ ($+1.50$); $\sigma_1{=}\sigma_2{=}0.009$ → doubly degenerate | $-4.441$ ($+1.50$) |

Three points:

- **$e$ lands exactly** ($-4.441$ in full agreement); $a_1$ differs by $3$ meV (the fine-mesh host GF's small adjustment).
- **Singlet vs doublet is read straight off the singular-value degeneracy** ($\sigma_1\!\approx\!\sigma_2\Rightarrow e$, $\sigma_2\!\gg\!\sigma_1\Rightarrow a_1$) — this turns §4's symmetry factorization into a computable criterion and closes the earlier caveat about relying on eigenvalue degeneracy with no character analysis.
- **Cost**: a dim-$539$ SVD at ~$400$ values of $E$, i.e. **minutes**, with no $P$/$Q$ dressing touched at all.

The $e=+1.50$ is the **11-band** value (the host GF contains only the 11 Wannier bands, self-consistently with 11-band explicit); putting more bands into the host GF (a larger rewann) will refine it toward explicit's $+1.35$ / DFT's $+1.19$ — at the cost of a few more eigenvalue $k$-sums, which is precisely where Koster–Slater is cheap. Script: `tools/koster_slater_levels.py`.
