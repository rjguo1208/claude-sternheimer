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

## 5. Head-to-head: the previous production method ($\chi$ augmentation) on the same system

The pre-Sternheimer production recipe — full diagonalization of the explicit
$\le$150-band $H_{\rm eff}$ with **multi-center $\chi$ augmentation** (the same
10-center set that produced the historical $+217\to+41$ meV SOC calibration
ladder: 3 Mo + 7 S cage, Gram-truncated) — rerun on the *identical* non-SOC
6$\times$6 relaxed V$_S$ inputs (same cubes, same $\omega_0$), via EDI-direct +
$\chi$ against a commensurate-grid twin of the same NSCF.

![VS chi vs two-level](../assets/vschi_levels.png)

| method | gap $e$-doublet (eV, rel. VBM) | error | basis engineering |
|---|---|---|---|
| truth | $+1.1726$ (2$\times$) | — | — |
| explicit-150, no $\chi$ | $+1.3516$ (2$\times$) | $+179$ meV | none (the phantom) |
| explicit-150 + $\chi$, keep=30 | $+1.2077$ (2$\times$) | $+35$ meV | 10-center $\chi$ set + keep |
| explicit-150 + $\chi$, keep=38 | $+1.2038$ (2$\times$) | $+31$ meV (saturated) | 10-center $\chi$ set + keep |
| **MODE B (two-level)** | $+1.183$ / $+1.212$ | $+10$ / $+39$ meV | **none** |

Read: on the vacancy class the two cures are **equally accurate** (both land the
doublet in the $\sim$30 meV static-$\omega_0$ band; centroids $+33$ vs $+25$ meV).
The differences are structural: $\chi$ preserves the $e$-degeneracy exactly but
needs a defect-specific orbital set and a truncation parameter scanned to
saturation; the two-level ladder needs *no defect-specific input at all* — and it
also fixes the VBM-edge window ($\chi$: uniform $\sim+30$ meV overshoot on all
four edge states; the no-$\chi$ a$_1$ phantom $+0.296\to+0.089$/$+0.091$ under
either cure vs truth $+0.060$). Non-SOC replication of the historical ladder:
$+179 \to +31$ meV here vs $+217 \to +41$ with SOC.

## 6. O$_S$ and Se$_S$: where the tail helps — and where it cannot

Same MODE B chain (relaxed geometry, ghost 4-species NSCF for the foreign-element
nonlocal channel), truth-gated. Se$_S$ (isovalent, host-like states):

![SeS levels](../assets/ses_levels66.png)

| Se$_S$ | VBM edge (eV) | CBM edge |
|---|---|---|
| truth | $-0.016$ (2$\times$) / $-0.011$ | $+1.660$ (2$\times$) |
| MODE B | $+0.010$ (2$\times$) / $+0.017$ | $+1.688$ (2$\times$) |

A uniform $+27$ meV common-mode; **internal spacings reproduced to $\le 2$ meV**,
level ordering (doublet below single) restored where M-only inverts it, gap clean.

O$_S$ (a new element — O-2p orbitals absent from the pristine Bloch basis):

![OS levels](../assets/os_levels66.png)

| O$_S$ | O-2p pair | a$_1$ | mid-gap | CBM edge |
|---|---|---|---|---|
| truth | $-0.002$ (2$\times$) | $+0.082$ | **empty** | $+1.653$ (2$\times$) |
| M-only | $-0.006$ (2$\times$) | $+0.048$ | $+0.801$ (2$\times$) | $+1.659$ (2$\times$) |
| MODE B | $-0.054$ (2$\times$) | — | $+0.807$ | $+1.609$ (2$\times$) |

The $\approx+0.80$ eV mid-gap feature has **no truth counterpart** and appears in
every truncated-active variant; rest dressing at frozen $\omega_0$ moves it by
only $+6$ meV. This is the known O-2p **missing-weight push-up** (tiny weight
$\times$ 100-eV tails) in downfolded form: the state is *rest-dominant*, so a
static-$\omega_0$ self-energy — no matter how exactly the rest is inverted —
cannot relocate its active-space shadow. The complementarity rule from the SOC
campaign is thus reproduced inside one framework: **vacancy/displacement defects
$\to$ the Sternheimer tail cures them with zero basis engineering; new-element
defects $\to$ they need the foreign atom's $\chi$ in the diagonalized space (or an
$\omega$-resolved $\Sigma$)**. The two methods are not competitors but the two
halves of one toolbox.

## 7. Residuals and next levers

Remaining deviations are small and mechanistically assigned: the V$_S$ doublet
splitting (29 meV where the truth is degenerate) and $\sim$40 meV valence-edge
offsets, plus the Sdisp resonance region at $\sim$30–50 meV — all consistent with
the frozen-$\omega_0$ (static) fold evaluated 1.2–1.9 eV away from the states in
question. The designed fix is an $\omega$-resolved rest (Lanczos continued
fraction over the same operators); second lever, raising the explicit boundary
(NSCF bands 150 $\to$ 300) which lowers both the ladder contraction and the tail
weight.
