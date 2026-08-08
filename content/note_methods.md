# Method derivations (compact)

All methods diagonalize the same active-space effective Hamiltonian and differ
only in the rest self-energy $\Sigma$:

$$H_{\rm eff}(\omega)=\varepsilon_A+\frac{M_{AA}+\Sigma(\omega)}{N_k},$$

with $\varepsilon$ the pristine NSCF eigenvalues, $M_{mn}=\langle\psi_m|\Delta
V|\psi_n\rangle$ the (supercell-summed) defect matrix elements, $N_k$ the
coarse-grid count. On the $N_k$-commensurate grid the eigenvalues compare
level-by-level with the defect supercell — the truth gate.

## 1. Exact downfolding identity

Split the Hilbert space $\mathcal{H}=A\oplus R$ ($A$ = the $n_{\rm sub}$-band
manifold per $k$; $R$ = everything else, semicore + high bands + the full
plane-wave tail). Feshbach/Löwdin partitioning of $(\omega-H)^{-1}$ gives the
exact eigenvalue condition

$$\det\!\big[\,\omega-H_{AA}-\Sigma(\omega)\,\big]=0,\qquad
\Sigma(\omega)=V_{AR}\,\big(\omega-H_{RR}\big)^{-1}V_{RA},$$

with $H_{RR}=P_R(H^0+\Delta V)P_R$. Every method below is one approximation to
$\Sigma$. Paths that re-enter $A$ mid-way are excluded from $\Sigma$ by
construction — the diagonalization of $H_{\rm eff}$ generates them (no double
counting).

## 2. M-only

$$\Sigma=0.$$

## 3. Bare second order (frozen $\omega_0$)

One bounce through the explicit rest bands, $\Delta V$ to second order:

$$\Sigma_{ab}=\frac{1}{N_k}\sum_{n\in R,\,n\le 150}
\frac{M_{an}M_{nb}}{\omega_0-\varepsilon_n}.$$

## 4. Two-level static (MODE A/B)

All orders in $\Delta V$, still frozen at $\omega_0$. Split
$R=R_1\oplus R_2$ ($R_1$: explicit bands $\le 150$; $R_2$: the physical tail)
and Schur-factorize $(\omega_0-H_{RR})^{-1}$:

$$\Sigma(\omega_0)=V_{A1}\,D_1\,V_{1A}
+\tilde s_b^{\dagger}\,G_{22}\,\tilde s_a,$$

$$D_1=(\omega_0-H_1-V_{11})^{-1}\ \text{(exact, zheevd)},\qquad
\tilde s_a=P_2V\,(1+D_1P_1V)\,|\psi_a\rangle,$$

$$G_{22}=\big[\omega_0-H_2^0-W_{22}\big]^{-1},\qquad
W_{22}=V_{22}+V_{21}D_1V_{12},$$

with $G_{22}$ applied by the Neumann ladder
$x^{(n+1)}=D_2\,[\tilde s+W_{22}\,x^{(n)}]$, $D_2$ a projected CG solve.
Convergence requires $\rho(D_2W_{22})<1$ **in every channel**: a dressed rest
state near $\omega_0$ (a $\rho\to1$ collective mode) makes a finite-rung ladder
silently incomplete — the MODE B caveat documented in the results erratum.
MODE A replaces the physical $R_2$ by explicit model bands (dense algebra,
exact reference gate).

## 5. $\omega$-resolved block Lanczos (MODE C)

No split, no ladder, no frozen frequency. Sources
$\chi_a=P_R\,\Delta V\,|\psi_a\rangle$; block-tridiagonalize the full dressed
rest operator with block width $N_A$:

$$H_{RR}Q_j=Q_{j-1}B_{j-1}^{\dagger}+Q_jA_j+Q_{j+1}B_j,\qquad
Q_0R_0=\mathrm{qr}(\chi),$$

then for **any** complex $\omega$ the resolvent follows from the block
continued fraction (evaluated bottom-up, $S_N=(\omega-A_N)^{-1}$):

$$S_j=\big[\omega-A_j-B_j^{\dagger}S_{j+1}B_j\big]^{-1},\qquad
\Sigma(\omega)=\frac{1}{N_k}\,R_0^{\dagger}\,S_1(\omega)\,R_0.$$

The inverse is constructed, never expanded — no convergence-radius condition;
new poles of the dressed $H_{RR}$ (collective defect modes) are carried
exactly. Levels are quasiparticle fixed points and spectra come free:

$$e_i(\omega)=\omega,\qquad
A(\omega)=-\tfrac{1}{\pi}\,\mathrm{Im}\,\mathrm{Tr}
\big[\omega+i\eta-H_{\rm eff}(\omega+i\eta)\big]^{-1}.$$

One chain serves every $\omega,\eta$; the static methods are the single-point
evaluation $\Sigma(\omega_0)$ of the same object.

## 6. $\chi$ augmentation (basis route)

The alternative philosophy: do not fold $R$ into $\Sigma$ — put the missing
directions into the diagonalized space. Bloch-summed atomic orbitals
$\chi_{lm}(k)$ at chosen centers, orthogonalized and Löwdin-whitened against
the bands, extend the basis:

$$H=\begin{pmatrix}\varepsilon+M/N_k & M_{b\chi}/N_k\\
M_{\chi b}/N_k & h'+M_{\chi\chi}/N_k\end{pmatrix},\qquad
h'=\langle\tilde\chi|H^0_{\rm pristine}|\tilde\chi\rangle,$$

diagonalized in full. Cures missing-weight push-ups by supplying the
directions themselves; requires a defect-specific orbital set (and a
truncation rank scanned to saturation).

## 7. Practical requirement: the fold grid

$M_{mn}=\langle\psi_m|\Delta V|\psi_n\rangle$ is evaluated by folding the
supercell $\Delta V$ onto the primitive grid,

$$V_q(\mathbf r)=\sum_{\mathbf R}\Delta V(\mathbf r+\mathbf R)\,
e^{i\mathbf q\cdot(\mathbf r+\mathbf R)},$$

implemented as a real-space modulo map $r^{\rm cube}_i \to r^{\rm cube}_i \bmod
n^{\rm prim}_i$. It is exact **only** if

$$n^{\rm cube}_i = N^{\rm sc}_i \; n^{\rm prim}_i \quad\text{(supercell multiplicity, not merely an integer ratio)} .$$

A divisor that is not the multiplicity fails just as badly: folding a 240-point
cube of a $6\times$ supercell onto a 30-point primitive grid (ratio 8, perfectly
integral) treats the cell as eight primitive cells and returns 15 spurious states
in the gap window where there should be 5, with every degeneracy broken.
Otherwise the map mis-assigns positions: $\Delta V$
is aliased, C$_3$ degeneracies split by tens of meV and $\lVert\tilde
V\rVert$ can be off by tens of percent — with in-gap states still appearing at
roughly the right energies, so the failure is silent. A 6$\times$6 cube of
$240\times240\times300$ therefore needs a primitive NSCF with an explicit
`nr1=nr2=40, nr3=300`. EDT prints the two grids and the detected multiplicity on the first fold and
**aborts** unless $n^{\rm cube}=N^{\rm sc}n^{\rm prim}$ holds on all three axes;
degeneracy splitting of a symmetric defect is the cheap independent check.

The pair may equally be satisfied from the other side, and more cheaply: band-limit
the cube onto $162\times162\times216$ by G-space truncation and pair it with the
natural ecut-50 grid $27\times27\times216$. Truncation is exact for every Fourier
component the coarse grid holds — here it discards $4\times10^{-8}$ of the weight —
while the fold cost falls as $n^{\rm cube}$: $5.6$ s per Lanczos step instead of
$15.6$, with quasiparticle levels moved by $\le 0.1$ meV.

### 7a. 2D slabs: skipping the vacuum in the contraction

Most of a slab cell is vacuum, and $V(\mathbf r)\psi(\mathbf r)$ there is
multiplying zeros. The FFTs still need the whole grid, but the real-space
contraction does not — and with one rank per pool the $z$-slices are contiguous
in the flat index, so restricting it is a block list, not an indexing rewrite.

What a dropped slice costs a matrix element is bounded by
$\sum_{\mathbf r\in\rm drop}|V(\mathbf r)|\rho(\mathbf r)$, so **that product**,
not the density alone, is the profile to threshold. The distinction is not
cosmetic: thresholding $\rho$ alone keeps $87\%$ of the grid at $10^{-5}$ and
saves nothing, because in a slab the upper conduction bands are vacuum
resonances and $\rho$ has a vacuum plateau. Weighted by $|V|$ (`zslab_tol`):

| `zslab_tol` | $z$-slices kept | operator unit test | fold (s/step) | level shift |
|---|---|---|---|---|
| off | 216 (100%) | $1.9\times10^{-14}$ | 5.6 | — |
| $10^{-5}$ | 181 (84%) | $5.2\times10^{-8}$ | 5.50 | $0.000$ meV |
| $10^{-4}$ | 158 (73%) | $2.6\times10^{-7}$ | 5.38 | $0.000$ meV |
| $\mathbf{10^{-3}}$ | **95 (44%)** | $4.8\times10^{-4}$ | **4.90** | $\mathbf{0.005}$ **meV** |
| $10^{-2}$ | 46 (21%) | $3.7\times10^{-2}$ | 4.50 | $6.3$ meV |
| $5\times10^{-2}$ | 30 (14%) | $1.0\times10^{-1}$ | 4.40 | $17$ meV |

$10^{-3}$ is the knee: $44\%$ of the grid for a $0.005$ meV shift. Note that the
Lanczos matrices themselves move much earlier ($\lVert\Delta A\rVert$ reaches
$0.4$ on a scale of $31$ at $10^{-3}$) while the levels do not — the basis
rotates but spans the same space, which is why the *levels*, not the matrices,
are the acceptance test.

The saving is real but capped, and the profile says why: after the coarse cube,
$V\!\cdot\!\psi$ is only $1.7$ s of the $5.6$ s fold, and $10^{-3}$ cuts it to
$0.9$ s — exactly proportional to the grid, but $12\%$ of the step. The largest
remaining item is the KB nonlocal projector work ($1.3$ s), which is where the
next factor has to come from.

### 7b. The $k$-grid is free: the supercell is a box, not a periodicity

The supercell exists to *capture one isolated defect*, not because we want a
defect array. EDI builds the fold accordingly — on the fly, at the **exact**
$q=\mathbf k_f-\mathbf k_i$, with no BZ wrap, and with the phase running over
the whole cube (`ed_coarse_full_q`):

$$V_q(\mathbf r)=\sum_{\mathbf R\in\rm cube}\Delta V(\mathbf r+\mathbf R)\,
e^{i\mathbf q\cdot(\mathbf r+\mathbf R)} ,$$

which is the finite-box transform of a *localized* $\Delta V$. It is defined,
and generically nonzero, for **arbitrary** $q$ — so $\langle\psi_{\mathbf
k}|\Delta V|\psi_{\mathbf k'}\rangle$ connects any pair of $k$-points and the
$k$-grid need not match the supercell. (Wrapping $q$ into the first BZ *without*
the compensating phase is what would be wrong; MODE C's canonical storage
applies it explicitly, $V_{q+G}(\mathbf r)=e^{iG\cdot\mathbf r}V_q(\mathbf r)$,
exact because $e^{iG\cdot\mathbf R}=1$.) Measured on a $12\times12$ grid against
a $6\times6$ cube: matrix elements between $k$-points differing by a
non-supercell vector reach $0.975$ of the largest within-supercell one. Treating
the cube as a periodicity would discard most of the scattering channels.

Only the real-space modulo map constrains the grids, and that constraint is
§7's $n^{\rm cube}=N^{\rm sc}n^{\rm prim}$.

### 7c. Two different $12\times12$ calculations

Because $H_{\rm eff}=\varepsilon_A+(M_{AA}+\Sigma)/N_k$ carries the
Born&ndash;von-K&aacute;rm&aacute;n normalization, $N_k$ *is* the dilution: one
defect per $N_k$ primitive cells. Two inequivalent things therefore go by
"$12\times12$":

| | $k$-set | defect density | supercell truth? |
|---|---|---|---|
| **coset-decomposed** | four shifted $36$-$k$ runs | $1/36$ cells | yes — the $6\times6$ cell, at four supercell $k$ |
| **direct** | one $144$-$k$ run | $1/144$ cells | no — would need a $431$-atom cell |

The coset route splits the grid into $(N_k/N^{\rm sc})^2$ shifted runs of the
designed size, each of which folds to one supercell $k$-point. It holds the
concentration fixed, which is exactly what makes it comparable to a supercell
calculation — it is a *validation* device, not the physics target. The direct
run is the physics: full $q$ resolution at a four-times more dilute defect, the
limit EDI is built for. The two agree where they must — the coset-diagonal
blocks of the direct $144$-$k$ vertex reproduce the four shifted $36$-$k$ runs
to $2\times10^{-6}$ on a block scale of $66$ (independent NSCF runs, so this is a
gauge-invariant comparison of block eigenvalues).

## 8. From the downfold to the $T$-matrix and the electron–defect self-energy

The downfold output is the **partially scattered vertex** — rest space summed to
all orders, active space not yet iterated:

$$\tilde V(\omega)=\frac{M_{AA}+\Sigma(\omega)}{N_k},\qquad
H_{\rm eff}(\omega)=\varepsilon_A+\tilde V(\omega).$$

With the free active propagator $G^0(\omega)=\big[\omega+i\eta-\varepsilon_A\big]^{-1}$
(diagonal), the active-space multiple scattering is resummed by the
Lippmann–Schwinger equation

$$T(\omega)=\tilde V(\omega)+\tilde V(\omega)\,G^0(\omega)\,T(\omega)
\;\;\Longrightarrow\;\;
\boxed{\,T(\omega)=\big[\,1-\tilde V(\omega)G^0(\omega)\,\big]^{-1}\tilde V(\omega)\,}$$

so that the dressed propagator and the spectral function are

$$G(\omega)=G^0+G^0TG^0=\big[\,\omega+i\eta-\varepsilon_A-\tilde V(\omega)\,\big]^{-1},
\qquad A(\omega)=-\frac{1}{\pi}\,\mathrm{Im}\,\mathrm{Tr}\,G(\omega).$$

The last identity is the bookkeeping check: $T$ is a repackaging of the same
$H_{\rm eff}$, so levels, DOS and $T$-matrix poles are one object. Division of
labour: rest scattering lives in $\Sigma$, active scattering in the $[1-\tilde
VG^0]^{-1}$ resummation — never both (no double counting).

**Electron–defect self-energy.** In the dilute single-site limit at defect
density $n_d$ (defects per primitive cell) the disorder-averaged self-energy of
a Bloch state is the diagonal $T$-matrix element,

$$\Sigma^{\rm ed}_{n\mathbf k}(\omega)=n_d\,T_{n\mathbf k,\,n\mathbf k}(\omega),$$

and its imaginary part gives the quasiparticle linewidth, equivalently (optical
theorem) a golden rule with $T$ in place of the bare matrix element:

$$\frac{1}{\tau_{n\mathbf k}}=-\frac{2}{\hbar}\,\mathrm{Im}\,
\Sigma^{\rm ed}_{n\mathbf k}(\varepsilon_{n\mathbf k})
=\frac{2\pi}{\hbar}\,n_d\!\!\sum_{m\mathbf k'}\!
\big|T_{n\mathbf k,\,m\mathbf k'}(\varepsilon_{n\mathbf k})\big|^{2}\,
\delta(\varepsilon_{n\mathbf k}-\varepsilon_{m\mathbf k'}),$$

with the momentum-relaxing (transport) rate obtained by inserting the usual
$\big(1-\mathbf v_{n\mathbf k}\!\cdot\!\mathbf v_{m\mathbf k'}/|\mathbf
v_{n\mathbf k}|^{2}\big)$ weight. Setting $T\to\tilde V\to M_{AA}/N_k$ recovers
the standard Born-level electron–defect matrix elements — the $T$-matrix is
their all-order completion.

Two practical notes. **(i)** As computed on an $N_k$ grid the construction
describes exactly one defect per $N_k$ cells, $n_d = 1/N_k$; a different dilute
density means re-solving $T$ on a finer grid (the vertex itself is
density-independent). **(ii)** $\tilde V(\omega)$ is smooth in $\mathbf k$
(defect-localized, hence short-ranged in $\mathbf R,\mathbf R'$) and is the
object to Wannier-interpolate; $T$ is *not* smooth — its bound-state poles and
resonances must be generated by solving $[1-\tilde VG^0]^{-1}$ on the fine grid.
With MODE C the vertex carries its true frequency dependence, so $T(\omega)$ has
the correct analytic structure; static vertices pinned at $\omega_0$ displace
features by tens of meV once $|\omega-\omega_0|\sim 1$ eV.

## 9. Does the vertex survive Wannier interpolation?

Section 7c leaves the method with a hard limit: the fold and the
reorthogonalization both scale as $N_k^3$ at fixed rank count, and $N_k$ *is* the
defect concentration, so a finer mesh is simultaneously more expensive and a
different physical system. The way out is the Lu&ndash;Bernardi construction —
carry the vertex in a **pair-Wannier** basis,

$$\mathcal M(\mathbf R_e,\mathbf R_p)=\frac1{N_k^2}\sum_{\mathbf k\mathbf k'}
e^{+i\mathbf k\cdot\mathbf R_e}e^{-i\mathbf k'\cdot\mathbf R_p}\,
U^\dagger(\mathbf k)\,M(\mathbf k,\mathbf k')\,U(\mathbf k') ,$$

and read it back at any $\mathbf k$. That is only legitimate if $\mathcal M$
decays in **both** arguments, so the question has to be measured, not assumed.

### 9a. Both indices are electron positions; the defect is the origin

$\mathcal M(\mathbf R_e,\mathbf R_p)=\langle w_{\mathbf R_e}|\Delta V|w_{\mathbf
R_p}\rangle$ — **both** labels are Wannier lattice vectors of the *electron*.
The defect is not a third index. It is where $\Delta V$ lives, hence the point
the kernel decays about, and in EDI's own supercells it sits at the origin so it
never appears explicitly. Ours sits at the cell centre, $\mathbf R_{\rm
def}=(3,3)$ in primitive lattice units, which is why it does. The conventional
way to remove it is a translation,

$$\Delta V(\mathbf r+\mathbf R_{\rm def})\;\Longleftrightarrow\;
M(\mathbf k,\mathbf k')\to e^{i(\mathbf k-\mathbf k')\cdot\mathbf R_{\rm def}}M(\mathbf k,\mathbf k') ,$$

after which the ordinary origin-centred Wigner&ndash;Seitz construction applies.
Choosing the WS representatives about $\mathbf R_{\rm def}$ instead, as done
here, is the same thing.

Three implementations agree on this, up to a change of variables. EDI's code
transforms the two indices separately — `edbloch2wane` over the electron $k$ at
fixed $q$, `edbloch2wanr` over $q$ — the electron&ndash;phonon parametrization,
and it already writes a `decay.M` file carrying exactly the envelope diagnostic
plotted below. The repository's `ft_convention.md` records the symmetric
Lu&ndash;Bernardi double transform, which is the form used here, and this
project's earlier `tmatrix_p6_wannier.py` uses that form with a defect-centred
truncation. The two parametrizations are a shear,
$\mathbf R_e^{\rm EDI}=\mathbf R_p-\mathbf R_e$ and $\mathbf R_p^{\rm EDI}=\mathbf R_p$:
EDI's first index is the *separation* of the two Wannier functions, ours is each
one's distance from the defect.

**Why the origin matters.** $\mathcal M$ is built on Born&ndash;von-K&aacute;rm&aacute;n
representatives, and images differ by $N^{\rm sc}\mathbf a$. On the coarse grid
$e^{i\mathbf k\cdot N^{\rm sc}\mathbf a}=1$, so the choice is invisible — the
transform inverts exactly whatever you pick. On a *finer* grid it is not $1$, and
a wrong representative aliases: measured, that turns a $1$ meV prediction into a
$300$ meV one. Note the asymmetry — truncating a coarse-grid kernel (what
`tmatrix_p6_wannier.py` does) is insensitive to this; only interpolation is.

$\mathbf R_{\rm def}$ is read from the structure — the S present in the pristine
cube and absent from the defect cube. Taking it from the peak of the kernel
diagonal also works, but returns a representative valid **only on the lattice it
was found on**: $(3,3)\equiv(-3,-3)$ modulo 6 but not modulo 12, so the same
"empirical centre" is right for a $6\times6$ R-lattice and wrong for a
$12\times12$ one. That is the mistake that produced the $300$ meV run here.

### 9b. Leave-one-out on the vertex

The direct $144$-$k$ run (section 9 of the results note) supplies $M(\mathbf
k,\mathbf k')$ on all $144^2$ pairs in one Wannier gauge — the same
wavefunctions are wannierized, `exclude_bands = 1-6, 18-20` leaving exactly the
11 active bands 7&ndash;17 as an isolated composite group, so $U$ is unitary
($\lVert U^\dagger U-1\rVert=2\times10^{-10}$, $\Omega_D=0.03$ &Aring;$^2$).
Build the kernel from the $\Gamma$ coset alone — a complete $6\times6$ dataset,
$6.2\%$ of the pairs — and predict the rest.

![Wannier interpolation of the defect vertex](../assets/wannier_decay.png)

Panel (a) is every one of the $2\,982\,529$ matrix elements against the radius at
which a spherical truncation would drop it; the envelope falls from $1.46$ on
site to $2\times10^{-4}$ at $6a$ and the median element is four decades below it
throughout. Panel (b) is the diagonal $\mathcal M_{mm}(\mathbf R,\mathbf R)$ of
each Wannier orbital, which shows what the defect actually couples to: the three
S $p$ orbitals on the vacancy sublattice start at $1.43$&ndash;$1.46$, the five
Mo $d$ at $0.05$&ndash;$0.10$, and the far S $p$ at $0.003$ — and all of them are
below $10^{-6}$ by $6a$.

| pairs | count | median | p90 | max |
|---|---|---|---|---|
| both on the $6\times6$ grid (training) | 1296 | $4\times10^{-15}$ | — | $1\times10^{-14}$ |
| one off-grid | 7776 | $0.64$ meV | $0.95$ | $1.53$ |
| **both off-grid** | **11664** | $\mathbf{0.74}$ **meV** | $1.03$ | $2.44$ |

(Element-wise, in the same gauge; quoted as $\delta M\,\mathrm{Ry}/N_k$, the
scale at which a vertex element enters $H_{\rm eff}$.) The training pairs come
back at machine precision, which certifies the transform rather than the
physics; the verdict is the $93.8\%$ that were never seen, and they land at
about **1 meV**.

Panel (a) says why. The kernel falls from $1.46$ on site to $2\times10^{-4}$ at
$6a$, and the $6\times6$ and $12\times12$ curves — the same quantity measured in
two different BvK cells — lie on top of each other wherever they overlap, so the
smaller cell is not aliasing, only truncating early. A $12\times12$ coarse grid
buys another factor of four of decay beyond what a $6\times6$ already delivers.

Set against the other errors in this project — a $+7.6$ meV downfold residual and
$4$&ndash;$11$ meV of periodic-image error in the $6\times6$ reference itself —
interpolation would not be the dominant error. What this does **not** yet test is
the $\omega$-dependent part: only the bare $M_{AA}$ was interpolated here.
$\Sigma^R(\omega)$ should be at least as local for $\omega$ in the gap, since it
carries an extra exponentially decaying resolvent whose range is set by the
distance from $\omega$ to the rest spectrum — but that is an argument, and the
same leave-one-out is the measurement.

## 10. Summary

| method | $\Sigma$ | order in $\Delta V$ | $\omega$ |
|---|---|---|---|
| M-only | $0$ | 0 | — |
| bare-150 | one rest bounce | 2 | frozen $\omega_0$ |
| MODE A/B | Schur two-level + ladder | all (per channel: rung-limited) | frozen $\omega_0$ |
| MODE C static | block-Lanczos CF at $\omega_0$ | all, exact | frozen $\omega_0$ |
| **MODE C $\omega$-res** | block-Lanczos CF | all, exact | **resolved** |
| $\chi$ aug. | none (basis enlarged) | all within the enlarged basis | — |

Numbers and gates: see the [Two-level results](twolevel-results.html) page.
