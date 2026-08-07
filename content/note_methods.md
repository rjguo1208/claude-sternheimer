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
n^{\rm prim}_i$. It is exact **only** if $n^{\rm cube}_i/n^{\rm prim}_i$ is an
integer in every direction. Otherwise the map mis-assigns positions: $\Delta V$
is aliased, C$_3$ degeneracies split by tens of meV and $\lVert\tilde
V\rVert$ can be off by tens of percent — with in-gap states still appearing at
roughly the right energies, so the failure is silent. A 6$\times$6 cube of
$240\times240\times300$ therefore needs a primitive NSCF with an explicit
`nr1=nr2=40, nr3=300`; the ecut-50 automatic grid ($27\times27\times216$) is
not usable. EDT prints the two grids on the first fold and warns when the ratio
is not integral; degeneracy splitting of a symmetric defect is the cheap
independent check.

## 8. Summary

| method | $\Sigma$ | order in $\Delta V$ | $\omega$ |
|---|---|---|---|
| M-only | $0$ | 0 | — |
| bare-150 | one rest bounce | 2 | frozen $\omega_0$ |
| MODE A/B | Schur two-level + ladder | all (per channel: rung-limited) | frozen $\omega_0$ |
| MODE C static | block-Lanczos CF at $\omega_0$ | all, exact | frozen $\omega_0$ |
| **MODE C $\omega$-res** | block-Lanczos CF | all, exact | **resolved** |
| $\chi$ aug. | none (basis enlarged) | all within the enlarged basis | — |

Numbers and gates: see the [Two-level results](twolevel-results.html) page.
