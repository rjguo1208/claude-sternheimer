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
invert only on the defect's compact support. **This requires $\tilde V^W(R',R)$ to be short-ranged.**

**The Wannierization itself is correct.** Reconstructing $H^W(R)$ from the same $U(k)$ pipeline
used for $\tilde V$ reproduces the independently-built `hr.dat` $H_W(R)$ to $2.3\times10^{-3}$ eV
(Fig. 1, left), and $H_W(R)$ decays cleanly ($\sim23$ eV at $R=0$ to $<0.01$ eV by $|R|=5$) — the
Wannier functions are well localized and the gauge is consistent.

![Left: gauge test — H_W(R) from hr.dat vs the P5-b pipeline (overlap, clean decay). Right: the downfolded-potential coupling from the origin cell, which does not decay.](../assets/vtilde_locality.png)

*Figure 1. **Left:** the gauge test — $H_W(R)$ from `hr.dat` (blue) and reconstructed from the
P5-b rotation (orange) coincide and decay $\sim3000\times$ over six cells, so the rotation/FT is
correct. **Right:** the coupling $\lVert O^W(R',R{=}0)\rVert$ for the Born $M$, rest $\Sigma$, and
$\tilde V$ — all essentially **flat** in $|R'|$, i.e. not localized.*

**But $\tilde V^W$ does not localize for this data.** $M^W$, $\Sigma^W$ and $\tilde V^W$ are flat
in $R$ (Fig. 1, right), and a truncated Koster–Slater inversion only converges to the full result
when the cutoff covers nearly the entire grid — no speed-up. The cause is visible in the bare
matrix element $M(k_f,k_i)$ itself (Fig. 2): it is $\approx f(k_f-k_i)$ (translation-invariant to
$77\%$) with weight concentrated on **commensurate momentum transfers** — the diagonal stripes are
the Bragg structure of a *periodic vacancy array*.

![Structure of the Born matrix element M(k_f,k_i): diagonal stripes (left), the k_f map at k_i=Gamma (middle), and |M| collapsed onto q=k_f-k_i (right).](../assets/vtilde_qstructure.png)

*Figure 2. The Born matrix element $\lVert M(k_f,k_i)\rVert$ shows diagonal stripes (left) — the
hallmark of an $M\!\approx\!f(q)$, periodic-potential-like operator — with $|M(q)|$ broad across
the BZ and peaked on commensurate transfers (right), not a single smooth $q\!=\!0$ form factor.*

**Diagnosis.** The $V_d/V_p$ cubes are a **$6\times6$ supercell** (cell $\approx36$ Bohr $\approx6\times$
the primitive $5.97$ Bohr). The vacancy potential does not fully decay within $6\times6$, so the
periodic defect images interact and $M$ carries the array's commensurate Bragg structure;
Wannierized on the $12\times12$ Born–von-Kármán grid, $\tilde V^W$ inherits the $6\times6$
periodicity and looks delocalized. This is **not** a long-range / range-separation problem (the
defect is neutral and $|M(q)|$ is broad, not the $1/q$ divergence of a charged defect) and **not**
a gauge bug (the gauge test passes). It is a **supercell-size** effect.

**Fix.** A larger defect supercell (so $\Delta V$ decays inside it $\Rightarrow$ smooth isolated-defect
$M(q)\Rightarrow$ localized $\tilde V^W$) is needed for the Koster–Slater speed-up; a denser $k$-grid
alone does not help (the $6\times6$ periodicity persists). The coarse-grid direct inversion (P5-a)
needs no locality and is unaffected.

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

P0–P3 and P5-a are complete and validated; P5-b is implemented and has produced a clear,
actionable finding. Open items:

- [ ] **P6 — transport:** feed $|T_{PP}(\omega)|^2$ on-shell into the golden-rule rate (replacing
  EDI's $|M|^2$) for a beyond-Born vs first-Born mobility on the coarse grid (uses P5-a directly).
- [ ] **larger defect supercell** to localize $\tilde V^W$ and unlock the Koster–Slater / fine-grid
  route (P5-b).
- [ ] **T9 / T5:** rest-BZ-grid convergence of $\tilde V$ and Wannier-gauge invariance.
