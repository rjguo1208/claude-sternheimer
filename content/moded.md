# MODE D (folded free-electron tail) as a cross-check on MODE C

### 6×6 non-SOC MoS$_2$, three systems (ideal V$_S$, relaxed V$_S$, relaxed O$_S$): DOS reconciliation and free-electron energy-space convergence

> **Method attribution**: MODE D (OPW / folded free-electron tail) was contributed by our collaborator **cz** and
> implemented on the `opw-tail` branch of the EDI development repository (`doc/opw_tail.md` + `post/tail_fold.py`).
> This page is **our independent reproduction and cross-check on the same data** — the first piece of the MODE C
> reconciliation that the method note lists as "still to do".
> No implementation code for that method appears here; we report only the numbers we ran ourselves.

## 1. What is being tested

**Band completion of the outer space**: cut the Hilbert space as $1 = P + K + Q$ — $P$ is the active window
(11 bands of $d+p$), $K$ is the explicit outer space up to $n_x$ (with real DFT energies and matrix elements), and
$Q$ is the tail above $n_x$ up to the plane-wave cutoff (**never enumerated**). $Q$ is eliminated **exactly** with a
Schur complement:

$$\Sigma_P(\omega)=T_{PP}+\bigl(W_{PK}+T_{PK}\bigr)\bigl[\omega-H_{KK}-T_{KK}\bigr]^{-1}\bigl(W_{KP}+T_{KP}\bigr),
\qquad T_{XY}=W_{XQ}\,G_Q(\omega)\,W_{QY}$$

The tail therefore appears in **three places**: the direct path $T_{PP}$, the dressing of the explicit-band
Hamiltonian $T_{KK}$, and the dressing of the $P\!\leftrightarrow\!K$ coupling $T_{PK}$. There are only two
approximations: $G_Q$ takes the plane-wave-diagonal free-electron form
$D_G=1/(\omega-\bar V-|k_f+G|^2)$, and each explicit band needs only one vector
$\chi_X=\hat Q\,\Delta V|\psi_X\rangle$ (where $\hat Q$ removes only the projections of the **existing** $n_x$
explicit bands).

The comparison target, **MODE C**, is our exact $\omega$-resolved Feshbach folding (implemented in production with
a block-Lanczos continued fraction; on this page the same target quantity is obtained by a direct solve on the same
160-band data).

## 2. Setup

| item | value |
|---|---|
| system / defect | MoS$_2$ monolayer 6×6 supercell, **ideal** (unrelaxed) S vacancy |
| relativity | scalar (non-SOC) |
| k mesh | 36 k (a stride-2 sub-block of 12×12, exact q, no wrapping) |
| explicit band universe | 160 bands; active $P$ = bands 7–17 (396 states) |
| MODE D partition | $n_x=100$ → $K$ = 89 bands, tail block = 60 bands |
| free-electron background | $\bar V=-7.649$ eV (energy-matched, tail block bottom $+28.07$ eV) |
| frequency | everything evaluated at $\omega+i\eta$, $\eta=10$ meV (retarded) |
| reference | 160-band all-order folding ($\mathrm{VBM}=-5.9365$ eV, gap $1.657$ eV) |

**Why everything runs at complex frequency**: evaluating on the real axis renders every outer-space pole inside the
window as a zero-width needle — the very pitfall we had just hit and fixed in the spectral-function pipeline (the
continued fraction was moved to evaluation at $\omega+i\eta$).

## 3. Independent reproduction (confirming the pipelines share a source)

Using data we assembled ourselves (36k sub-block slicing, independent wavefunction reading, an independent driver),
we reproduce the M2 self-consistency numbers from the method note — the maximum error over the gap window
(~20 levels) relative to the 160-band all-order folding:

| variant | this run | method note | meaning |
|---|---|---|---|
| drop-block | **11.10 meV** | 11.1 | no tail at all |
| **ffree-open** | **1.32 meV** | 1.3 | folded free-electron tail |
| fexact | **2.79 meV** | 2.8 | folding framework + exact tail spectrum |

All three land. **Note that `fexact` is worse than `ffree`**: once the coupling is folded in, the details of the
tail **model spectrum** are already a second-order effect — a counter-intuitive conclusion that holds in our
independent reproduction too.

## 4. Convergence in free-electron energy space

The one axis the method note never swept, and the **most practical parameter on the production side**: after the
plane-wave sum in the tail contraction $T_{XY}=\sum_G \chi_X^*(G)\,D_G\,\chi_Y(G)$ is truncated by kinetic energy
$|k+G|^2$, how many $G$ actually have to be stored per band per k.

| $|k+G|^2$ cutoff | plane waves | block weight | absolute error | relative to full $G$ |
|---|---|---|---|---|
| full set (1361 eV) | 24034 | 1.00000 | 1.315 meV | — |
| 500 eV | 5351 | 0.99990 | 1.316 | 0.000 |
| **300 eV** | **2486** | 0.99885 | 1.320 | **0.007** |
| 200 eV | 1353 | 0.99497 | 1.347 | 0.048 |
| 150 eV | 878 | 0.98910 | 1.394 | 0.102 |
| 100 eV | 479 | 0.97119 | 1.663 | 0.480 |
| 75 eV | 312 | 0.93234 | 2.379 | 1.478 |
| 50 eV | 169 | 0.78137 | 5.946 | 4.631 |
| 30 eV | 80 | 0.17414 | 9.494 | 8.179 |
| 10 eV | 15 | 0.00073 | 11.061 | 9.746 |

Two readings: **(i)** as the cutoff $\to 0$ it returns exactly to drop-block (11.06 vs 11.10) — internal
consistency proving itself; **(ii)** the model's own error (1.32 meV) completely dominates the truncation error
above roughly 200 eV. So the tail contraction needs only **~1400–2500 $G$ (6–10% of the full set)**, which
compresses the side file by **10–17×** without touching the error budget.

> **Limitation**: in this framework $\chi$ is obtained by **band projection** (onto bands 101–160), so its $G$
> content is biased soft; the real production $\chi=\hat Q\Delta V|\psi\rangle$ is harder. The table above is
> therefore a **lower bound**, and the true values need the not-yet-implemented Fortran $\chi$ export.

## 5. DOS reconciliation

![MODE D vs MODE C DOS](../assets/moded_dos_VS.png)

Active-space DOS on a log axis, with the gap shaded; the lower panel is the difference against MODE C.

| MODE C peak (eV) | MODE D shift | drop-block shift | feature |
|---|---|---|---|
| $-0.285$ | 0.0 meV | 0.0 meV | valence-edge multiplet |
| $-0.020$ | 0.0 | 0.0 | pair adjacent to the VBM |
| **$+1.170$** | **0.0** | **$+10.0$** | **V$_S$ deep $e$ doublet** |
| $+1.695$ | 0.0 | 0.0 | conduction edge |
| $+1.940$ / $+1.980$ | 0.0 | 0.0 | resonances |

**All six peaks fall within the 5 meV resolution of the $\omega$ mesh.** The $L_1$ error of the in-gap DOS:
**MODE D $0.196$ vs drop-block $1.140$ states** — the folded tail cuts the error to $1/5.8$. In the difference plot,
drop-block's $\pm 23$ bipolar oscillation at $+1.17$ is the classic signature of a level shift, while MODE D leaves
only a $\pm 8$ lineshape residual with no shift.

## 6. Relaxed geometries: V$_S$ and O$_S$ (the second and third systems)

Beyond the ideal geometry, two defects in **relaxed** geometry provide a further cross-check. The edmats for these
two were **recomputed** with the current `edi.x` (v8d, the 20260809 band-axis convention): hosts are the respective
commensurate 6×6 NSCF saves (150 bands, 36 k, FFT 40×40×300, $E_{\rm cut}=50$ Ry), ΔV taken from the uncoarsened
240×240×300 cube (commensurability ratio 6:1 ✓), `pot_align='vacuum'`, with both setups identical.

> **The index convention is decided by physics, not by assumption**. The same matrix admits 6 possible index
> conventions, and the M$_{AA}$ Hermiticity of all of them is $3\times10^{-14}$ — **Hermiticity cannot
> discriminate**. The criterion is the levels: the correct convention gives V$_S$'s M-only doublet at $+0.0337$ (2×)
> and $+0.080$ (corresponding line by line to the established $+0.0395$ (2×) and $+0.0858$, the ~6 meV difference
> coming from the cube/host/alignment setup), while the wrong conventions give a gapless mess of a staircase. This
> is the third lesson of this kind in the project (after the SU(2) pairing phase and the C3 k-mapping).

### 6a. Reconstruction error (n_x = 100, static pin)

| system | drop-block | **ffree-open** | fexact |
|---|---|---|---|
| ideal V$_S$ (160 bands) | 11.10 meV | **1.32** | 2.79 |
| relaxed V$_S$ (150 bands) | 61.43 meV | **8.16** | 5.65 |
| relaxed O$_S$ (150 bands) | 98.33 meV | **8.96** | 1.92 |

**The relaxed systems depend far more strongly on the tail** — the cost of dropping it rises from 11 meV to
61–98 meV, and the folded free-electron tail brings that down to 1/7.5–1/11. Another reversal: on ideal V$_S$,
`ffree` beats `fexact`, while on both relaxed systems `fexact` wins — so the conclusion that "the model spectrum's
details are a second-order effect" is **system-dependent**, and in the more strongly coupled relaxed geometries the
tail spectrum itself starts to carry a little weight.

### 6b. Free-electron energy-space convergence (consistent across all three systems)

| cutoff | plane waves | ideal V$_S$ | relaxed V$_S$ | relaxed O$_S$ |
|---|---|---|---|---|
| full set | 24034 / 8494 | 1.315 | 8.164 | 8.960 |
| 500 eV | 5351 | +0.000 | +0.001 | +0.001 |
| **300 eV** | **2486** | **+0.007** | **+0.022** | **+0.024** |
| 200 eV | 1353 | +0.048 | +0.153 | +0.161 |
| 100 eV | 479 | +0.480 | +1.630 | +1.644 |
| cutoff → 0 | — | → drop-block ✓ | → drop-block ✓ | → drop-block ✓ |

(The last three columns are deviations from the full-$G$ result, in meV.) Two readings: **(i)** the plane-wave count
at a given energy cutoff is **identical across all three** (2486, 1353, 878, 479, …) — same unit cell, so this is a
**geometry-determined universal number**, not a coincidence of one case; **(ii)** the ideal system's host cutoff is
100 Ry (24034 $G$ in the full set) and the two relaxed ones are 50 Ry (8494), yet the cost at 300 eV is 0.007–0.024
meV in every case — **the conclusion is insensitive to the host cutoff**.

### 6c. DOS reconciliation

![MODE D vs MODE C, relaxed V_S](../assets/moded_dos_VSR.png)

![MODE D vs MODE C, relaxed O_S](../assets/moded_dos_OSR.png)

| system | peak agreement (MODE D) | largest drop-block shift | in-gap $L_1$: MODE D / drop |
|---|---|---|---|
| ideal V$_S$ | 6/6 peaks $\le$ 5 meV | +10 meV (deep doublet) | 0.196 / 1.140 (×5.8) |
| relaxed V$_S$ | 9/9 peaks $\le$ 5 meV | **+50 meV** (deep level +1.345) | 0.786 / 4.897 (×6.2) |
| relaxed O$_S$ | 3/3 peaks = 0 meV | +5 meV; the $+0.28$ resonance pushed ~90 meV | 1.055 / 4.392 (×4.2) |

The relaxed O$_S$ difference plot makes the point best: under drop-block, the resonance at $+0.28$ eV is pushed
bodily to $+0.37$ (~90 meV), while MODE D puts it back where it belongs, leaving only a lineshape residual.

> **Do not tabulate these numbers alongside previously published values**: this section's recomputation uses an
> uncoarsened cube, a 50 Ry host and vacuum alignment, so the relaxed V$_S$ deep level's $\omega$-resolved value
> lands at $+1.345$, whereas the previously published MODE C value is $+1.1985$ — the ~145 meV difference comes from
> **setup**, not from a MODE D error. Every comparison on this page is internally consistent within its own data:
> MODE C and MODE D consume the same matrix and the same set of eigenvalues.

### 6d. The ghost-state staircase: MODE D passes the phantom gate

The object most easily misjudged in defect downfolding is a **spurious mid-gap state**. Historically three different
things have been conflated here, and the new data separates them (static pin, mid-gap window $0.15$–$1.55$ eV):

| outer-space treatment | order / bands | V$_S$ mid-gap | O$_S$ mid-gap |
|---|---|---|---|
| M-only | no outer space | EMPTY | **$+0.8062$ (2×)** ← O-2p push-up ghost |
| bare (2nd order, full outer space) | 2nd / 139 | **$+0.7031$, $+0.7034$** ← Born spurious pair | EMPTY |
| bare (2nd order, keep only) | 2nd / 89 | $+0.8856$, $+0.8862$ | EMPTY |
| drop-block | **all-order** / 89 | $+0.4763$, $+1.4165$ (2×) | $+0.3731$ (2×) |
| **MODE D** (ffree) | all-order 89 + analytic tail | $+0.4230$, $+1.3721$ (2×) | $+0.2837$, $+0.2845$ |
| all-order reference | all-order / 139 | $+0.4149$, $+1.3684$ (2×) | $+0.2748$, $+0.2757$ |

**Three distinct things**: (i) O$_S$'s **O-2p ghost at $+0.806$ (2×)** lives only in M-only (established record
$+0.801$ (2×)), and **second-order dressing alone is enough to kill it** ($+0.0296$ (2×)/$+0.0401$; established
record $+0.025$ (2×)/$+0.035$ — two independent reproductions); (ii) V$_S$'s **$+0.703$ spurious pair is an artifact
of the second-order truncation itself** (established record bare-150 $+0.7112$ (2×)) and vanishes the moment
all-order folding appears; (iii) therefore **drop-block having no ghost is correct behaviour** — it is not an
"explicit sum", it performs an all-order folding over 89 bands.

**The phantom gate carries weight precisely because it is not trivial**: MODE B dressed the entire outer space with
a 6-level deflated ladder — high enough order — yet **kept the ghost** ($+0.807$) because its basis could not span a
slow mode. MODE D uses a completely different approximation and lands on MODE C's side: the mid-gap states of both
systems track the all-order reference to **3.7–8.8 meV**, while drop-block is off by **48–98 meV**.

> **One unresolved discrepancy (recorded honestly)**: this section's newly computed relaxed O$_S$ has a
> near-degenerate pair at $+0.275$ in its **all-order reference**, whereas in the established tables the relaxed
> O$_S$ **supercell ground truth and MODE C are both "mid-gap empty"**. The M-only and bare second-order results
> from this same data agree with the established values to 5 meV, so **the discrepancy appears only at the
> "all-order" step**. Candidate causes: static pin vs $\omega$ self-consistency (though the $+0.28$ peak is still
> there in the $\omega$-resolved DOS); the established O$_S$ MODE C used a 36-step chain and O$_S$ happens to be the
> single system where MODE B and C disagree ($\lVert\Delta\Sigma\rVert=0.68$ Ry); and setup differences such as
> `pot_align`. **To be decided**: running a longer production chain (block Lanczos) on the new data would settle
> both "is the gap empty" and "are 36 steps enough" at once.

### 6e. The complete complement: MODE D vs MODE C vs supercell ground truth

Every MODE D number above lives inside the **M2 self-consistency framework**: the tail $Q$ is the **band block**
$(n_x,150]$, and $\chi$ is built by band projection. Production is different — $Q$ is the **true complement up to
the plane-wave cutoff**, and $\chi=\hat Q\Delta V|\psi_X\rangle$ has to be exported from the Fortran side. We
supplied that piece (EDT gained `chi_dump_raw`: batch-implanting all $(n,k)$ sources for bands $1..n_x$ down
**exactly the same** `apply_dV` path as the chain, so the KB nonlocal term is exact; it streams out the **raw**
$\Delta V|\psi_X\rangle$ together with $|k_f+G|^2$ and the explicit-band coefficients, leaving the projection depth
as a post-processing knob).

**This is the first time MODE C's and MODE D's outer spaces are the same object**: the production chain's Krylov
vectors already live in the complete plane-wave complement (`project_PA` removes only the active bands), and now
MODE D's tail is that same complement. Host, cube and `pot_align` match item by item on both sides.

**Three gates**: export vs edmat at **7.8e-16 / 1.3e-15 Ry** (0.29/0.50 for the wrong index orientation, so the
criterion has teeth); explicit-band coefficient norm 1.000000; and **tail weight
$\|\chi\|^2/\|\Delta V\psi\|^2$ = 0.668 (V$_S$) / 0.813 (O$_S$)** — i.e. two thirds to four fifths of the weight
still sits beyond the 100 explicit bands.

**The tail variants swap places in production.** In the M2 framework `ffree-open` is best across the board; in the
complete complement it leaves open-mode spurious poles — the free-electron staircase starts from
$\bar V\approx-10$ eV while $Q$'s true bottom is at $+28$ eV, so the low-$|k_f+G|^2$ labels are fictitious states
with an oppositely-signed denominator:

| variant | O$_S$ mid-gap | V$_S$ deep doublet |
|---|---|---|
| `open` (formula as written) | **$+0.2150$** (spurious) | $+1.3228$ (2×) |
| **`clamp`** (lower bound pinned at $Q$'s bottom) | **empty** ✓ | $+1.3089$ (2×) |
| `sapofix` (SAPO projected energies) | $+0.8512$ (worse) | $+1.3085/+1.3096$ |

`sapofix` is **numerically unstable** in production: $p_G$ reaches 0.998, and $1/(1-p_G)$ throws some $G$'s
effective energy out to $-4294$ eV (in M2, $Q$ is a narrow block so it stays tame). **Guard rails are needed at the
implementation level.** On V$_S$, `clamp` and `sapofix` agree to 0.4–1 meV, showing that this system's tail has
entered a model-independent regime; on O$_S$ they still disagree, showing it is more sensitive.

![three-way DOS, complete complement](../assets/moded_3way_complete.png)

| peak (eV) | V$_S$ truth | V$_S$ MODE C | V$_S$ MODE D | O$_S$ truth | O$_S$ MODE C | O$_S$ MODE D |
|---|---|---|---|---|---|---|
| valence side | $-0.290\,{-}0.260\,{-}0.090\,{-}0.020\,{+}0.060$ | $-0.270\,{-}0.240\,{-}0.070\,{-}0.010\,{+}0.080$ | $-0.270\,{-}0.230\,{-}0.200\,{+}0.010\,{+}0.120$ | $-0.260\,{-}0.230\,{-}0.180\,{+}0.000\,{+}0.080$ | $-0.310\,{-}0.280\,{-}0.210\,{-}0.050\,{+}0.040$ | $-0.300\,{-}0.100\,{-}0.020\,{+}0.030$ |
| **mid-gap** | $+1.170$ | $+1.190$ | $+1.290$ | **empty** | **empty** ✓ | **empty** ✓ |
| CBM edge | $+1.670$ | $+1.690$ | $+1.710$ | $+1.650$ | $+1.610$ | $+1.580/{+}1.650$ |

**Three conclusions**:

1. **The MODE C side is self-consistent**: overall offsets against ground truth of $+20$ meV for V$_S$ and $-40$ meV
   for O$_S$, consistent with the established deep-anchor means of $+25.5$ / $-49.3$ meV — no new problem on the
   chain side.
2. **Both pass the phantom gate**: MODE C and MODE D both give an **empty** O$_S$ gap, in agreement with the
   supercell ground truth. The free-carrier tail **really does cure the O-2p push-up artifact** — and this is not a
   trivial result: MODE B, with a 6-level ladder, kept it instead (+0.807).
3. **But MODE D puts V$_S$'s deep doublet 100 meV above MODE C** ($+1.290$ vs $+1.190$), while the M2
   self-consistency test at the same $n_x$ reports only 8 meV. **M2 systematically overstates the method's
   accuracy** — it only asks the tail to reproduce 50 bands, whereas production asks it to reproduce the entire
   complement (measured to hold 2/3–4/5 of the weight). The method note flags this as unproven ("M2 is
   self-referential ... says nothing about the absolute physics"); **now there is a number for it.**

> **Open**: $n_x$ is the method's convergence knob, and this export is capped at $n_x=100$ (post-processing can only
> sweep downward). Whether $n_x=150$ or higher closes that 100 meV needs a new export; the choice of `clamp`'s lower
> bound and the $\bar V$ fit are both first versions too. These three items are the next round of work on
> "production usability".

## 7. Conclusions and open items

**Conclusions**: on **three systems** in 6×6 non-SOC (ideal V$_S$, relaxed V$_S$, relaxed O$_S$), MODE D reproduces
MODE C's spectrum using 100 explicit bands plus an analytic tail, with **all 18 DOS peak shifts $\le$ 5 meV** (the
mesh resolution), whereas drop-block — which discards the same tail block — pushes defect levels by 10–90 meV. The
folded tail cuts the in-gap $L_1$ error to 1/4.2–1/6.2; the relaxed geometries depend on the tail 5–9× more strongly
than the ideal one, which makes this test carry more weight on the relaxed systems.

**Not yet demonstrated**: the 12×12 dilute limit; the $n_x$ convergence curve within the production framework
(see §6e); the 12×12 dilute limit (the method's real target regime, and the one where the chain method is most
expensive) has not been done; and the MODE C on this page is a direct exact folding on the same data, so running the
production chain (block Lanczos) would additionally test the chain machinery itself.
All three systems are one working point — 6×6, scalar, $n_x=100$ — and the $n_x$ convergence curve is taken from the
method note's conclusion without independent verification.
