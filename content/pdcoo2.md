# PdCoO$_2$ delafossite: band structure (PBE and PBE+U)

### A standalone DFT calculation starting from one POSCAR: two cells, two high-symmetry paths, and the $U_{\mathrm{Co}\text{-}3d}$ sensitivity

> This page is independent of the rest of the site (defect $T$-matrix downfolding); it is a separate band-structure
> calculation. Inputs and scripts are in `~/edi_tmatrix/pdcoo2/`, and the runs were done on a single Anvil highmem node.

## 1. Identifying the structure

The given POSCAR is an **R$\bar 3$m rhombohedral primitive cell** (4 atoms = 1 formula unit):

$$a_{\rm rho}=6.1359\ \text{Å},\qquad \alpha=26.666^\circ$$

Converting via $a_{\rm hex}=|\mathbf a_1-\mathbf a_2|$ and $c_{\rm hex}=|\mathbf a_1+\mathbf a_2+\mathbf a_3|$,

$$a_{\rm hex}=2.8300\ \text{Å},\qquad c_{\rm hex}=17.7430\ \text{Å}$$

i.e. the standard cell of **PdCoO$_2$ delafossite** — the quasi-two-dimensional oxide famous for having a
**room-temperature in-plane resistivity lower than metallic copper**.
Pd sits at 3a, Co at 3b (CoO$_6$ octahedra, Co$^{3+}$ low-spin $d^6$, hence **non-magnetic**), and O at 6c ($z=0.1112$).

**The hexagonal conventional cell is constructed numerically from the rhombohedral one** (no Wyckoff table is
applied): take $\mathbf a_1-\mathbf a_2$, $\mathbf a_2-\mathbf a_3$ and $\mathbf a_1+\mathbf a_2+\mathbf a_3$ as
the basis, fill with rhombohedral translations and deduplicate, giving 12 atoms; the result lands automatically on
the R-centring translations $(0,0,0)$, $(2/3,1/3,1/3)$, $(1/3,2/3,2/3)$ — **the construction proves its own
correctness through the result**.

## 2. Computational setup

| item | value | rationale |
|---|---|---|
| pseudopotentials | PseudoDojo NC-SR v0.5 PBE stringent (Pd 18e / Co 17e / O 6e) | the same set as everything else on this site |
| cutoff | `ecutwfc` = 100 Ry | the order suggested for Co-3d stringent |
| occupations | `smearing`, Marzari–Vanderbilt, `degauss` = 0.02 Ry | **metal** |
| spin | `nspin` = 1 | Co$^{3+}$ low-spin $d^6$; non-magnetic in both experiment and calculation |
| $+U$ | `HUBBARD (ortho-atomic)`, $U_{\mathrm{Co}\text{-}3d}=4$ eV | a common value; compared point by point against PBE |
| SCF $k$ mesh | rhombohedral $12^3$; hexagonal $16\times16\times3$ | the hexagonal cell is scaled anisotropically for $c/a=6.27$ |
| band path | hexagonal $\Gamma$-M-K-$\Gamma$-A-L-H-A (281 points); rhombohedral RHL1 (301 points) | RHL1 parameters computed from the actual $\alpha$: $\eta=0.8206$, $\nu=0.3397$ |

## 3. Bands

![PdCoO2 bands](../assets/pdcoo2_bands.png)

Top: hexagonal conventional cell (12 atoms); bottom: rhombohedral primitive cell (4 atoms). Solid blue is PBE,
dashed red is PBE+U, and $E_F$ is aligned to zero.

## 4. Three readings

**(i) Only one band crosses the Fermi level.** The rhombohedral primitive cell counts **1**; the hexagonal cell
counts 3 — exactly the same band folded into a 3× cell (a self-consistency check). That single, wide, nearly free
electron-like Pd band is the origin of PdCoO$_2$'s extreme conductivity, and the Fermi surface is a simple
hexagonal cylinder.

**(ii) The quasi-two-dimensionality is immediately visible.** Along $\Gamma\to$A (out of plane) the band is almost
**perfectly flat**, while the in-plane dispersion along $\Gamma$-M-K-$\Gamma$ reaches $\sim 8$ eV — this is the
band-structure explanation for the three-orders-of-magnitude anisotropy in resistivity.

**(iii) $+U$ acts in a very clean place, and the transport physics is insensitive to it.**

| quantity | PBE | PBE+U (4 eV) | difference |
|---|---|---|---|
| $E_F$ (hexagonal) | 14.8272 eV | 14.8304 eV | **3 meV** |
| the Pd band near the Fermi level | — | — | essentially coincident |
| Co-3d occupied bands ($-1.5$ to $-4$ eV) | — | — | pushed down by $\sim 1$ eV overall and rearranged |

The reason is structural: the states near the Fermi level are **Pd-4d/5s** in character, while $U$ acts only on
**Co-3d**; and the $t_{2g}$ of Co$^{3+}$ low-spin $d^6$ is a **filled** shell (closed shell), so the Hubbard term
only shifts a filled band rigidly and cannot change the topology near the Fermi level.
**Conclusion: for transport and the Fermi surface, PBE suffices; $+U$ is needed only to place the Co-3d occupied
states (for comparison with XPS/XAS).**

## 5. Two self-checks

| check | result |
|---|---|
| total energy per formula unit for the two cells (PBE) | $-645.26000$ (hexagonal, $-1935.78001/3$) vs $-645.26077$ (rhombohedral) → **0.77 mRy = 10 meV** |
| same ($+U$) | $-645.10444$ vs $-645.10407$ → **0.37 mRy** |
| number of bands crossing $E_F$ | rhombohedral 1 ↔ hexagonal 3 (3 f.u. folding) ✓ |

The two cells come from **completely independent** constructions (the rhombohedral one taken directly from the
POSCAR, the hexagonal one generated programmatically), with different $k$ meshes and different symmetry
identification, and their total energies agree to 10 meV/f.u.; with $+U$ they tighten to 0.37 mRy — which shows
that the Hubbard projection acts on the same set of Co-3d orbitals in both cells, so **the $+U$ implementation is
self-consistent too**.

## 7. Should $U$ be applied at all: the answer from linear response

The $+U$ in §4 used the common value of 4 eV. To turn "should it be applied" from a judgement call into a citable
computed quantity, we **computed $U$** with **DFPT linear response** (QE's `hp.x`, the Cococcioni–de Gironcoli
definition):

$$U=\left(\chi_0^{-1}-\chi^{-1}\right)_{\rm Co,Co}$$

| quantity | value |
|---|---|
| **$U_{\mathrm{Co}\text{-}3d}$ (linear response)** | **7.95 eV** |
| bare response $\chi_0^{\rm CoCo}$ | $-0.8022$ eV$^{-1}$ |
| self-consistent response $\chi^{\rm CoCo}$ | $-0.1085$ eV$^{-1}$ |
| screening ratio $\chi/\chi_0$ | **0.135** |

(ortho-atomic projection, $n_q=2\times2\times2$, `conv_thr_chi` = 1e-6, 8 virtual supercell points; 21 minutes
wall clock.) The self-consistent response is only **13.5%** of the bare one — the Co-3d occupations are strongly
pinned under perturbation, so **this is a genuinely localized shell**, and the computed $U$ is **twice** the
commonly used 4 eV.

### 7a. Magnetic ground-state check

"Closed shell" underpins every argument on this page, so it should be **demonstrated rather than assumed**. A
spin-polarized SCF (`nspin` = 2, Co starting moment 0.4) gives:

| quantity | value |
|---|---|
| total / absolute magnetization | $0.00$ / $0.00$ $\mu_B$/cell |
| total energy of the spin-polarized solution | $-645.26075840$ Ry |
| total energy of the non-magnetic solution | $-645.26076524$ Ry |
| difference | $7\times10^{-8}$ Ry $\approx 0.1$ μeV |

Starting from a finite moment, self-consistency **collapses exactly back to the non-magnetic solution**. The
low-spin $d^6$ closed shell $t_{2g}^6e_g^0$ is established.

### 7b. Overall judgement: $U$ is large, but should not be used for transport

Three facts have to be read together:

1. $U_{\rm Co}$ = 7.95 eV — Co-3d **really is localized**;
2. the moment is strictly 0 — the correlated shell is **closed**, and DFT+U's
   $\frac{U}{2}\mathrm{Tr}[\mathbf n(1-\mathbf n)]$ **vanishes identically** for an idempotent occupation matrix,
   leaving only the rigid shift of occupied states down and empty states up;
3. at $U$ = 4 eV, $E_F$ moves only 3 meV and the Pd band at the Fermi level barely moves (§4).

Raising $U$ from 4 to 7.95 eV only pushes the Co-3d occupied bands down about twice as far (~2 eV instead of
~1 eV) and **still leaves the region near the Fermi level untouched**. Therefore:

> **For transport, the Fermi surface, and electron–phonon: use plain PBE.**
> **To place the Co-3d occupied states (for XPS/XAS/RIXS comparison): use $U=7.95$ eV — the self-consistently
> computed value, not a guessed 4 eV.**
> For anything more rigorous, consider U+V (`hp.x` also gives the interlayer Co–O $V$) or HSE/GW.

**Not converged**: the $U$ above is at $n_q=2\times2\times2$. HP's $U$ has a q-mesh convergence, so a
$3\times3\times3$ point (about an hour) should be added before quoting it.

> **A practical trap** (worth recording): `hp.x` requires **all Hubbard atoms to come first in
> `ATOMIC_POSITIONS`**. We wrote them in POSCAR order as Pd-Co-O, so `hp_init` exited under control — but the
> actual message was completely buried under 48 ranks' worth of `forrtl: severe (28)` stack traces and could not be
> grepped out of the log. It was **the stack itself** (`hp_init.f90:50 → hp_stop_smoothly`) that identified this as
> a controlled exit rather than a crash, and reading those few lines of source found the cause. Reordering to
> Co-Pd-O passed on the first try.

## 9. Wannierization: 16 MLWFs reproduce the entire $p$-$d$ manifold

### 9a. Why 16, and why no disentanglement is needed

First the manifold structure of the bands (rhombohedral primitive cell, 301-point path, relative to $E_F$):

| bands | range (eV) | assignment |
|---|---|---|
| 9–10 | $-20.8 \ldots -18.9$ | O-2s |
| **11–26** | $\mathbf{-8.55 \ldots +2.28}$ | **O-2p(6) + Co-3d(5) + Pd-4d(5) = 16** |
| 24 | $-0.757 \ldots +0.952$ | the only conduction band crossing $E_F$ |
| 27+ | from $+3.23$ | free-electron-like high bands |

The gap between bands 10↔11 is **10.3 eV** and between 26↔27 is **0.95 eV** (cross-checked on a uniform mesh as
17.16→18.11 eV) — **the $p$-$d$ manifold is completely isolated**. So we take `num_bands = num_wann = 16` with
`exclude_bands = 1-10, 27-40` and **need no disentanglement at all**, which is the numerically most stable regime.

The orbital counting is self-consistent too: 16 bands = 6+5+5, and the two above $E_F$ (25, 26) are exactly the
Co-$e_g$ — the empty $e_g$ of Co$^{3+}$ low-spin $d^6$, corroborating the closed-shell conclusion of §7.

### 9b. Two sets of projections: $\Omega_I$ agrees to 9 digits

The `nscf` uses $8\times8\times8=512$ $k$-points (`nosym`/`noinv`) with plain PBE (see §7b). The two variants
differ only in the Pd projection:

| | projection | $\Omega_I$ | $\Omega_D$ | $\Omega_{OD}$ | $\Omega$ (Å$^2$) |
|---|---|---|---|---|---|
| **A** | Pd: $d$(5) + Co: $d$(5) + O: $p$(6) | 12.433543233 | 0.00999622 | 0.62411529 | **13.067655** |
| **B** | Pd: $s$(1) + Pd: $d_{xy},d_{xz},d_{yz},d_{x^2-y^2}$(4) + Co: $d$(5) + O: $p$(6) | 12.433543233 | 0.00999573 | 0.62410262 | **13.067642** |

$\Omega_I$ is a **gauge invariant**, and the two agree to **9 significant figures**; the total spreads differ by
$1.3\times10^{-5}$ Å$^2$. This is exactly how an isolated manifold should behave — the projection sets the starting
point of the optimization, not its endpoint. Swapping Pd's $d_{z^2}$ for an $s$ still converges to **the same global
minimum**, showing this MLWF solution is **not** the product of some local minimum.

The 16 WFs have a mean spread of **0.817 Å$^2$**, $\Omega_{OD}/\Omega=4.8\%$, and centres landing **exactly on the
atomic sites with zero drift**:

| WF | centre $z$ (Å) | assignment | spread (Å$^2$) |
|---|---|---|---|
| 1 | 0 | Pd $d_{z^2}$ | **1.528** |
| 2–5 | 0 | Pd $d_{xy},d_{xz},d_{yz},d_{x^2-y^2}$ | 0.755, 0.755, 0.948, 0.928 |
| 6–8 | 8.8715 $(=c/2)$ | Co $t_{2g}$ | 0.546, 0.540, 0.540 |
| 9–10 | 8.8715 | Co $e_g$ | 0.640, 0.640 |
| 11–13 | 1.98 $(=0.1112\,c)$ | O(1) $2p$ | 1.107, 0.758, 0.759 |
| 14–16 | 15.76 $(=0.8888\,c)$ | O(2) $2p$ | 1.107, 0.758, 0.759 |

Two pieces of physics that were **not designed in, but produced by the data itself**:

1. **Pd $d_{z^2}$ has a spread of 1.528 Å$^2$, twice that of the other four Pd-$d$.** It is precisely the orbital
   that makes up the conduction band — delocalized, strongly hybridized with Pd-5s, and therefore the least
   localized. This is a direct quantitative fingerprint of the $s$-$d$ hybridization picture of PdCoO$_2$'s
   conduction band.
2. **Co's five $d$ orbitals split automatically into 3 + 2** (three at 0.54, two at 0.64), exactly the octahedral
   crystal-field $t_{2g}$/$e_g$ splitting — even though the projections were not written out by symmetry.

### 9c. Point-by-point comparison against the DFT bands

![PdCoO2 Wannier bands](../assets/pdcoo2_w90_bands.png)

**The comparison method itself needs explaining**: we did **not** use wannier90's own `*_band.dat`. Its path
sampling does not necessarily coincide with the 301 points of the DFT band calculation, and the interpolation error
introduced by aligning the two is the same size as the quantity being measured. Instead we read `<seed>_hr.dat`
directly and built
$H(\mathbf k)=\sum_{\mathbf R}e^{2\pi i\mathbf k\cdot\mathbf R}H(\mathbf R)/\deg(\mathbf R)$
(551 Wigner–Seitz $\mathbf R$ vectors) **at the DFT's own 301 $k$-points**, then diagonalized — **rigorously
comparable point by point**.

| | RMS over all 16 bands | RMS for $|E-E_F|<1$ eV | max deviation | conduction band (24) RMS | conduction band max |
|---|---|---|---|---|---|
| **A** | **3.574 meV** | **2.825 meV** | 59.15 meV | **2.726 meV** | 13.31 meV |
| **B** | 3.574 meV | 2.825 meV | 59.15 meV | 2.724 meV | 13.31 meV |

On an $8^3$ mesh, the interpolation error over the whole 16-band manifold is **3.6 meV RMS**, **2.8 meV** near the
Fermi level, and **2.7 meV** for the conduction band (13 meV maximum). The thin lines in panel (c) are the
point-by-point errors of all 16 bands, mostly below 1 meV; the deviations concentrate near band crossings and can be
pushed lower by refining the $k$ mesh.

**An incidental self-consistency check**: this rerun's scf gives $E_F=14.8809$ eV, **exactly the same** as the
independent calculation in §3.

### 9d. Two methodological lessons

**(i) Band extrema must be measured on a uniform mesh, never along a high-symmetry path.**
We also tried a single-band model retaining only the conduction band ($\text{num\_wann}=1$), with the frozen window
set from the path readings "band 23 top $=-0.175$, band 25 bottom $=+0.891$" — and wannier90 reported at the 165th
$k$-point that the frozen window contained 2 bands while the target had only 1. Checking the true values on the
512-point uniform mesh: band 25's actual minimum is $+0.803$ eV, **88 meV lower than the path reading** — it simply
does not lie on any high-symmetry line. A high-symmetry path is a **set of measure zero** in the BZ, so reading
extrema from it is invalid in principle.
(The single-band model was subsequently dropped as unnecessary and is not part of this page's conclusions, but the
lesson stands.)

**(ii) Windows must be anchored to absolute band positions, not to $E_F$.**
The $8^3$ nscf gives $E_F=14.8310$ eV, **50 meV below** the $12^3$ scf's 14.8809 (mesh and smearing differences).
Writing a disentanglement window as "$E_F\pm$" makes it drift bodily by those 50 meV, while the entire purpose of
the window is to bracket **one specific band** with margins of only tens of meV. A constraint must be anchored to
the invariant it is meant to reproduce — here the relative band positions, not the Fermi energy.

## 10. Cost and reproducibility

Four SCF runs plus four band calculations (two cells × PBE/+U) take **under 6 minutes in total** on a single Anvil
highmem node (hexagonal SCF 80 s, rhombohedral SCF 45 s, bands 30–90 s each). Input cards, SLURM scripts and
plotting scripts are in `~/edi_tmatrix/pdcoo2/`; the run directories are
`/anvil/scratch/x-rg47749/pdcoo2h` (hexagonal) and `.../pdcoo2r` (rhombohedral).

The linear-response $U$ calculation lives in `/anvil/scratch/x-rg47749/pdcoo2u` (`hp_scf.in` / `hp.in` /
`mag_scf.in`, output `pdcoo2u.Hubbard_parameters.dat`).

The wannierization lives in `/anvil/scratch/x-rg47749/pdw90` (inputs and SLURM scripts mirrored in
`~/edi_tmatrix/pdcoo2/w90/`): scf 44 s + nscf (512 $k$) 52 s + about 3 minutes per variant (`pw2wannier90.x`
dominating), on a single node with 48 ranks × 2 threads. Comparison/plot script: `w90bands.py`.

> **The submission header must be complete**: `--ntasks=48 --cpus-per-task=2 --exclusive`, all three of them.
> With only `-n 48`, 96 threads run on 48 cores and the node is shared as well, taking the scf from
> **2.1 s/iteration to 50 s/iteration** — with **nothing unusual in the log at all**. The only way to find it is
> to compare `PWSCF ... WALL` against a known-good job.

**Directions to continue**: use this 16-band MLWF set for the **Fermi surface** (direct visualization of the
hexagonal cylinder) and for transport; refine the $k$ mesh to push the interpolation error into the sub-meV range;
projected density of states (to verify quantitatively that the Fermi surface is Pd's); a $U$ convergence point at
$n_q=3\times3\times3$; or add spin–orbit coupling.
