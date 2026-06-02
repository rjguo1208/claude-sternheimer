# Numerical results: downfolded potential and active $T$-matrix

These are the first end-to-end numerical results of the EDT package on the **MoS₂ sulfur
vacancy** (primitive NSCF, 150 bands, $12\times12$ $k$-grid; 11-band valence active manifold;
defect potential from the $V_d/V_p$ cubes, vacuum-aligned). They cover three milestones:
the full downfolded-potential block $\tilde V=M+\Sigma$ (**P3**), the active-space $T$-matrix
$T_{PP}=[1-\tilde V G^A]^{-1}\tilde V$ (**P5-a**), and the Wannier / Koster–Slater locality
study (**P5-b**). All matrix elements are validated against EDI's first-Born kernels; no raw
data is published — only the summary tables and figures below.

Reference energy $\omega_0=-1.0941$ eV (valence-band maximum), Sternheimer broadening
$\eta=0.05$ eV, active window $[-11.74,-1.04]$ eV, $N_k=144$.

## 1. The downfolded potential block $\tilde V$ (P3)

The rest-space (distant-band) dressing is folded into an effective active-space potential
$$\tilde V_{(m k_f),(n k_i)} = \underbrace{\langle\psi_{m k_f}|\Delta V|\psi_{n k_i}\rangle}_{M\ (\text{Born})}
 \;+\; \underbrace{-\tfrac{1}{N_k}\sum_{k'}\langle s_{m k_f}^{Q}|\chi_{n k_i}\rangle}_{\Sigma\ (\text{rest, Sternheimer})},$$
assembled over the **full BZ** with the $1/N_k$ measure. The full $1584\times1584$ block
($11$ bands $\times\,144$ $k$) was computed by a pool-parallel $k'$-sum (1 wholenode node,
36 MPI ranks = 36 pools, 2 h 11 m), each rank owning the full $G$-grid so `h_psi` runs natively
at its pool-local channel $k'$.

**Validation (all pass):**

| check | result | what it confirms |
|---|---|---|
| Born limit (T1) | source ket $=$ EDI $M$ to $2.3\times10^{-13}$ Ry | $\Delta V$-as-ket (local + KB nonlocal) exact |
| per-rank $H_0$ gate | $\langle\psi|H_0|\psi\rangle=\varepsilon$ to $6.0\times10^{-10}$ eV (all ranks) | Sternheimer operator setup correct |
| Hermiticity | $\lVert\tilde V-\tilde V^\dagger\rVert=9\times10^{-12}$ (pre-symmetrization) | $M,\Sigma$ independently Hermitian $\Rightarrow$ assembly correct |
| in-situ cross-check | band-17/$\Gamma$: $\tilde V_{nn}=-0.11712$ Ry $=$ single-rank diagonal to 6 digits | cross-pool MPI assembly correct |

**Physics — strongly beyond-Born.** At the VBM ($\Gamma$, top valence band):
$M_{nn}=+0.701$, $\Sigma_{nn}=-0.819$, $\tilde V_{nn}=-0.117$ Ry — the rest dressing is larger
than, and opposite in sign to, the bare Born coupling. Over the whole manifold:

| quantity | value | reading |
|---|---|---|
| $\lVert M\rVert_F,\ \lVert\Sigma\rVert_F,\ \lVert\tilde V\rVert_F$ | $367,\ 422,\ 100$ Ry | $\lVert\Sigma\rVert>\lVert M\rVert$; $\tilde V$ is $\sim\tfrac14$ of $M$ |
| $\lVert\Sigma\rVert/\lVert M\rVert$ | $1.15$ | rest dressing exceeds first Born |
| states with $|\Sigma_{nn}|>|M_{nn}|$ | $976/1584$ $(62\%)$ | majority strongly dressed (sign flips) |
| $\lVert\tilde V_{\rm offdiag}\rVert/\lVert\tilde V_{\rm diag}\rVert$ | $24.6$ | scattering is overwhelmingly $k$-changing |

The rest dressing **screens** the bare defect potential: the downfolded $\tilde V$ is $\sim3.7\times$
weaker in norm than $M$. First-Born (EDI) therefore overestimates the coupling here — exactly the
regime that motivates the $T$-matrix.

## 2. Active-space $T$-matrix $T_{PP}$ (P5-a)

With $\tilde V$ in hand, the active multiple scattering is resummed by one small inversion on the
coarse grid,
$$T_{PP}(\omega)=[\,1-\tilde V\,G^A(\omega)\,]^{-1}\tilde V,\qquad
  G^A_a(\omega)=\frac{1}{N_k}\,\frac{1}{\omega-\varepsilon_a+i\eta},$$
where the BZ measure $1/N_k$ lives in $G^A$ (the same $1/N_k$ that appears in $\Sigma$ — see the
[$k'$-normalization note](note-kprime-normalization.html)).

**Validation (code correct):** as the coupling is scaled to zero ($\lambda\to0$) or
$\eta\to\infty$ (so $G^A\to0$), the solver returns $T\to\tilde V$ to $\sim10^{-4}$.

**Result — beyond-Born cuts the resummed scattering.** Evaluated at the VBM:

| input | resummed $\lVert T\rVert_F$ | amplification |
|---|---|---|
| bare Born $M$ | $\lVert T_M\rVert=584$ Ry | $\times1.6$ vs $M$ |
| downfolded $\tilde V$ | $\lVert T_{PP}\rVert=297$ Ry | $\times3.0$ vs $\tilde V$ |

$T_{PP}=297$ Ry sits a factor $\sim2$ **below** the Born-input $T_M=584$ Ry, and the two differ by
$53\%$ — the rest dressing materially changes the scattering operator. A scan over $\omega$ shows
$\lVert T_{PP}\rVert/\lVert\tilde V\rVert\approx1$ across most of the window but **jumps to $4.4$ at
the VBM** — the active multiple scattering is *resonant at the band edge, where the carriers live*.

## 3. Wannier representation and locality (P5-b)

The textbook way to avoid the large $(N_b N_k)$ inversion on fine grids is the Koster–Slater /
defect-Green's-function trick: rotate $\tilde V$ to a localized Wannier basis,
$\tilde V^W(R',R)=\!\big[U^\dagger(k_f)\,\tilde V(k_f,k_i)\,U(k_i)\big]$ Fourier-transformed, and
invert only on the defect's compact support. This requires $\tilde V^W(R',R)$ to be short-ranged —
which in turn requires the Wannier rotation $U(k)$ to be in the **same Bloch gauge** as the
$\psi_{nk}$ that build $M$ (exactly the consistency EDI's `edbloch2wane` enforces by construction).

**A gauge mismatch — found and fixed.** Our first attempt reused the existing `filukk`, which had
been produced by a *separate* 17-band Wannier90 run, whereas `edt.x` computes $\tilde V$ from the
**150-band NSCF** evc. Different runs carry different per-$k$ Bloch phases (and degenerate-subspace
bases), so $U^\dagger M U$ was not smooth. The fingerprint: the electron-index decay of
$M^W(R_e;q)$ collapsed at $q=0$ (phase-insensitive) but stayed **flat for $q\neq0$**
(phase-sensitive) — a per-$k$ gauge mismatch, *not* a real non-locality. (Band/eigenvalue checks
pass either way, since eigenvalues are gauge-free, which is why this hid initially.)

The fix is to **re-Wannierize on the 150-band NSCF with the identical Wannier space** — same
projections (Mo:$d$, S:$p$), same windows, the same 11 valence bands 7–17 (`exclude_bands = 1-6,
18-150`) — giving a `filukk` consistent with the evc that build $M$ (its Wannier interpolation
reproduces the NSCF bands to $2\times10^{-5}$ eV).

![Electron-index decay of M before (left, q!=0 flat) and after (right, all q decay) the gauge fix.](../assets/vtilde_gauge_fix.png)

*Figure 1. Electron-index decay $\lVert M^W(R_e;q)\rVert$. **Left (old `filukk`, 17-band run):**
$q\!\neq\!0$ is flat — the gauge mismatch. **Right (new `filukk`, re-Wannierized on the 150-band
NSCF):** every $q$ now decays together by $\sim10^{3}\times$ over $\sim5$ cells.*

**With the consistent gauge, the Wannier / Koster–Slater route works.** $\tilde V^W$ is now
localized — peaked on the defect cell and decaying $\sim250\times$ — and the truncated inversion
$T=[1-\tilde V^W G^A]^{-1}\tilde V^W$ converges quickly:

| $R_{\rm cut}$ | subspace dim | $\lVert T\rVert$ old (mismatched gauge) | $\lVert T\rVert$ new (consistent gauge) |
|---|---|---|---|
| 2 | 275 | 0.057 | 0.057 |
| 3 | 539 | 0.19 | 1.72 |
| 4 | 891 | 0.59 | **2.06 (99.6%)** |
| 6 | 1584 | 2.06 | 2.06 |

![Koster-Slater truncation: ||T(Rcut)|| vs cutoff, converges by Rcut=4 with the consistent gauge.](../assets/vtilde_ks_converge.png)

*Figure 2. Koster–Slater truncation $\lVert T(R_{\rm cut})\rVert$ vs the cutoff radius (subspace
dimension under each tick). With the gauge-consistent `filukk_150` (blue) the inversion converges
by $R_{\rm cut}=4$ (dim 891, $\sim56\%$ of the full 1584); with the mismatched gauge (red) it only
reaches the full value at the full subspace. The localized $\tilde V^W$ is what makes the truncation
effective — the payoff of the gauge fix.*

So the neutral S-vacancy potential **is** short-ranged, as expected — the earlier "flat
$\tilde V^W$" was entirely the gauge mismatch, *not* a supercell-size or range-separation problem.
On this coarse $12\times12$ grid the converged subspace ($R_{\rm cut}\!=\!4$) is $\sim56\%$ of the
full one; on a finer $k$-grid the fixed $\sim$3–4-cell defect extent becomes a small fraction —
which is where the Koster–Slater speed-up pays off. **P3 (the $\tilde V$ block) and P5-a ($T_{PP}$
in the Bloch basis) use no Wannier rotation and were unaffected throughout.**

## 4. On the magnitude of $\tilde V$

The Frobenius norm $\lVert\tilde V\rVert_F=100$ Ry can look alarming, but the *individual* matrix
elements are small (RMS $\approx0.06$ Ry, on-site $\tilde V_{nn}\approx-0.12$ Ry): the norm is large
only because it root-sum-squares $\sim2.5$ million matrix elements. The raw eigenvalues
($\tilde V\in[-80.6,+26.5]$ Ry, $M$ up to $211$ Ry) are the discrete-basis values; the physical
scattering operator carries the BZ measure $1/N_k$ (it enters every $T$-matrix product through
$G^A$), so the physical eigenvalues are $\sim\lambda/N_k\approx0.5$–$1.5$ Ry — bounded by
$\max|\Delta V|\approx12$ Ry and giving a multiple-scattering parameter $\lVert\tilde V G^A\rVert\sim\mathcal O(1)$
(strong but sensible; the $T$-matrix resummation is essential, possibly resonant). See the
[$k'$-normalization note](note-kprime-normalization.html).

## 5. Status and next steps

P0–P3 and P5-a/b are complete and validated (the P5-b gauge consistency is now resolved with a
`filukk` re-Wannierized on the 150-band NSCF). Open items:

- [ ] **P6 — transport:** feed $|T_{PP}(\omega)|^2$ on-shell into the golden-rule rate (replacing
  EDI's $|M|^2$) for a beyond-Born vs first-Born mobility (uses P5-a directly).
- [ ] **fine-grid Koster–Slater** transport using the gauge-consistent `filukk` (P5-b machinery),
  where the localized $\tilde V^W$ delivers the inversion-size speed-up.
- [ ] **T9 / T5:** rest-BZ-grid convergence of $\tilde V$ and Wannier-gauge invariance.
