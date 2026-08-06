# Two-level dressing: results (Sdisp probe & relaxed V$_S$, 6$\times$6)

## 1. What is compared

All methods diagonalize the same effective Hamiltonian
$H_{\rm eff}=\mathrm{diag}(\varepsilon_A)+(M+\Sigma)\,R_y/N_k$ on the 6$\times$6 grid
(the grid whose $H_{\rm eff}$ *is* the supercell-array problem, so its eigenvalues
compare level-by-level with the supercell $\Gamma$ spectrum — the like-for-like
"truth gate"). The methods differ only in the rest-space self-energy $\Sigma$:

| method | $\Sigma$ |
|---|---|
| M-only | none (pure 11-band explicit block) |
| bare 2nd order | $V_{AR}\,G^0_R\,V_{RA}$ over the explicit $\le$150 bands |
| MODE A | $R_1$ ($\le$150) dressed to all orders in $\Delta V$; no genuine $>$150 tail |
| **MODE B (two-level full)** | $R_1$ exact ($D_1$, zheevd) **+ physical $>$150 tail via the Sternheimer ladder** (Schur-partitioned; see the *Deflated ladder* page) |

Setup: non-SOC, MoS$_2$ 6$\times$6, $\omega_0$ at the VBM, active = the 11-band
$d$+$p$ manifold (396 states), $R_1$ = the other 139 NSCF bands (5004 states),
ladder with 4–6 rungs (measured contraction 0.06–0.13 per rung).

## 2. S-displacement probe (weak, sharp perturbation)

One S displaced $+0.30$ Å along $z$. The exact gap is clean; the CBM doublet sits
at $+6.5$ meV (2$\times$).

![Sdisp level spectra](../assets/sdisp_levels66.png)

| method | CBM doublet (meV) | gap |
|---|---|---|
| truth | $+6.5$ / $+6.5$ | clean |
| M-only | $+19.3$ / $+19.4$ | stray state $+73$ meV |
| bare 2nd order | $+12.7$ / $+12.9$ | valence states pushed $+69$ meV into the gap |
| MODE A (no tail) | $+16.8$ / $+17.0$ | **spurious deep state at $-199$ meV** |
| **MODE B** | **$+8.2$ / $+8.5$** | **clean** |

## 3. Relaxed V$_S$ (strong defect — the original nemesis)

The sulfur vacancy with the relaxed C$_3$ shell, the system whose gap $e$-doublet
historically read $+217$ meV too high with a truncated basis and $+41$ meV with the
best multi-center $\chi$ augmentation.

![relaxed VS level spectra](../assets/vsrlx_levels.png)

| method | gap $e$-doublet (eV, rel. VBM) | error | gap |
|---|---|---|---|
| truth | $+1.1726$ (2$\times$) | — | clean |
| M-only | *missing entirely* | $>$0.5 eV | — |
| bare 2nd order | $\sim$1.16 among spurious states at $+0.20\ldots+1.42$ | — | flooded |
| **MODE B** | **$+1.183$ / $+1.212$** | **$+10$ / $+39$ meV** | **clean** |

No augmentation orbitals, no tunable hyperparameters: the $+10/+39$ meV accuracy
comes from first principles (exact $R_1$ inverse + the dressed physical tail),
matching or beating the best hand-tuned $\chi$ result on the defect class that
motivated the whole program.

## 4. Engineering numbers

| item | value |
|---|---|
| MODE A gate vs independent python Schur | $6.3\times10^{-16}$ Ry |
| ladder vs exact zgesv (model tail) | $2\times10^{-11}$ Ry |
| measured rung contraction (physical sources) | 0.06–0.13 (worst-mode $\rho$: 0.15–0.32) |
| full V$_S$ chain (SCF $\to$ cubes $\to$ block $\to$ MODE B $\to$ gate) | 40 min (1 node) |
| MODE B memory (v2: sliced channels + streamed $x$) | $\sim$5 GB/rank (was $\sim$30) |

## 5. Residuals and next levers

Remaining deviations are small and mechanistically assigned: the V$_S$ doublet
splitting (29 meV where the truth is degenerate) and $\sim$40 meV valence-edge
offsets, plus the Sdisp resonance region at $\sim$30–50 meV — all consistent with
the frozen-$\omega_0$ (static) fold evaluated 1.2–1.9 eV away from the states in
question. The designed fix is an $\omega$-resolved rest (Lanczos continued
fraction over the same operators); second lever, raising the explicit boundary
(NSCF bands 150 $\to$ 300) which lowers both the ladder contraction and the tail
weight.
