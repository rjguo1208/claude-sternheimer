# Active/Rest Partitioning of the Defect $T$-matrix
### cRPA-style downfolding, the coupling-second-order effective potential, and a $k$-decoupled Sternheimer solution

---

## 0. Overview and the one-paragraph summary

We compute the electron–defect $T$-matrix by splitting the host Green's function into an **active** block $A$ (bands near the band of interest / $E_F$, where dynamical effects are strong) and a **rest** block $R$ (distant bands, treated statically). This is the constrained-RPA (cRPA) logic transported from screened interactions to multiple scattering. The construction has two layers, and the central message of this note is that **they are cleanly separable and individually controlled**:

- **Layer 1 (rest):** the rest bands renormalize the bare defect potential $V$ into an effective, energy-dependent potential $\tilde V(\omega)$. Downfolding to the active space makes $\tilde V$ **exactly second order in the active–rest coupling** (Feshbach), so "second order" here is not an approximation — it is exact. The only genuine approximations at this layer are (i) optionally dropping the rest self-dressing $V_{QQ}$, and (ii) the static limit $\omega\to\omega_0$. We obtain $\tilde V$ by solving a **Sternheimer equation** whose operator is the periodic host Hamiltonian, hence **block-diagonal in crystal momentum** ($N_k$ independent small solves) — this is what avoids the $N_k\times N_k$ dense coupled system. Full rest dressing is recovered, *while staying $k$-decoupled*, by an iterative ladder (Neumann/Krylov) whose fixed point is exactly the dressed Sternheimer solution.
- **Layer 2 (active):** the full dynamical multiple scattering is resummed **exactly** by a single small matrix inversion in the active space, $T_{PP}(\omega)=[1-\tilde V\,G^A(\omega)]^{-1}\tilde V$.

Neither layer ever requires forming or inverting a dense $(N_kN_b)\times(N_kN_b)$ matrix.

> **Scope & status relative to the EDI code (read §16 first if you are implementing).**
> The production EDI code (`/anvil/projects/x-che190065/rjguo/qe-7.5/edi-dev`) currently computes the electron–defect matrix element strictly in the **first Born approximation** — a single defect scattering event, $M_{nk,mk'}=\langle\psi_{nk}|\Delta V|\psi_{mk'}\rangle$, fed into Fermi's golden rule. In the language of this note that is just the leading term $T\!\approx\!V$ (one power of the difference potential, no resolvent applied). **Nothing in this note is implemented yet:** there is no Sternheimer solve, no $G_0$ resolvent, no $T$-matrix resummation, and no $A/R$ projector anywhere in `src/`. The method here is the *beyond-Born* generalization of what EDI already does. The good news, detailed in **§16**, is that EDI's existing machinery lines up almost one-to-one with the abstract objects below — the Wannier active manifold *is* the active space $A$, the excluded bands *are* the rest $R$, the difference potential `V_colin`(+KB) *is* $V$, and the interpolated Wannier Hamiltonian *is* $G^A$ — so the construction can be built largely by reusing parts already present.

---

## 1. Setup and conventions

We work in the KKR / multiple-scattering convention. Let $H_0$ be the host (defect-free) crystal Hamiltonian, $G_0(\omega)=(\omega-H_0)^{-1}$ its retarded propagator, and $V=H-H_0$ the defect potential (in practice the supercell **difference potential**, local + nonlocal parts). The Lippmann–Schwinger equations are

$$
T(\omega)=V+V\,G_0(\omega)\,T(\omega),
\qquad
G(\omega)=G_0(\omega)+G_0(\omega)\,T(\omega)\,G_0(\omega).
$$

Formally $T=[1-VG_0]^{-1}V=V[1-G_0V]^{-1}$.

**Notation.** Throughout, $G^{A},G^{R}$ denote the relevant projections of the *host* propagator $G_0$ (not of $T$); the subscript $0$ is suppressed. We write $\omega_0$ for the reference energy (center of the active window, typically the band of interest or $E_F$), and $\Delta$ for the active–rest energy gap (the smallest $|\,\varepsilon_n-\omega_0|$ over $n\in R$).

| Symbol | Meaning |
|---|---|
| $P,\;Q$ | projectors onto active ($A$) and rest ($R$) subspaces, $P+Q=\mathbb 1$, $[P,H_0]=0$ |
| $G^A,\;G^R$ | active/rest projections of the host propagator $G_0$ |
| $H_0^Q=QH_0Q$ | rest-block host Hamiltonian |
| $V_{PP},V_{PQ},V_{QP},V_{QQ}$ | blocks $PVP$, $PVQ$, $QVP$, $QVQ$ of the defect potential |
| $\mathcal G^Q(\omega)$ | $V_{QQ}$-dressed rest propagator, Eq. (3.2) |
| $T^R$ | partial scattered $T$-matrix (rest-only resummation), §3 |
| $\tilde V(\omega)=PT^RP$ | active-space effective potential (downfolded), §3 |
| $T_{PP}$ | active block of the full $T$-matrix, §6 |

---

## 2. Active/rest partitioning (the cRPA structure)

Introduce orthogonal projectors $P$ (active) and $Q$ (rest), $P+Q=\mathbb 1$. Because the partition is by host band index and $[P,H_0]=0$, the host propagator is block-diagonal and separates additively:

$$
G_0(\omega)=G^A(\omega)+G^R(\omega),\qquad
G^A=\sum_{n\in A}\frac{|n\rangle\langle n|}{\omega-\varepsilon_n},
\quad
G^R=\sum_{n\in R}\frac{|n\rangle\langle n|}{\omega-\varepsilon_n}.
$$

$G^A=PG^AP$ lives entirely in the active space; $G^R=QG^RQ$ entirely in the rest space. The physical input — *strong dynamics in $A$, weak/static dynamics in $R$* — mirrors cRPA, where low-energy "target" screening is kept explicit and high-energy screening is folded into a partially screened interaction.

---

## 3. Partial scattered $T$-matrix (rest-only resummation)

Resum multiple scattering using **only** the rest propagator $G^R$:

$$
T^R(\omega)=V+V\,G^R(\omega)\,T^R(\omega)
\quad\Longrightarrow\quad
T^R=\big[1-VG^R\big]^{-1}V=V\big[1-G^RV\big]^{-1}.
\tag{3.1}
$$

This is the multiple-scattering analogue of the cRPA **partially screened interaction** $W^{r}=(1-vP^{R})^{-1}v$, under the dictionary $G^R\leftrightarrow P^{R}$, $V\leftrightarrow v$. It is the defect potential "pre-dressed" by the distant bands.

### 3.1 Downfolding to the active space = Feshbach effective potential

The object that enters the active dynamics is the active block $\tilde V\equiv PT^RP$. Expanding (3.1) and keeping the $PP$ block,

$$
\boxed{\;\tilde V(\omega)=PT^RP=V_{PP}+V_{PQ}\,\mathcal G^Q(\omega)\,V_{QP}\;}
\qquad
\mathcal G^Q(\omega)=\big[\omega-H_0^Q-V_{QQ}\big]^{-1}.
\tag{3.2}
$$

This is **exact** (Feshbach downfolding). The derivation is a direct resummation: in $T^R=V+VG^RV+VG^RVG^RV+\cdots$, every internal propagator is $G^R$, which lives in $Q$. Therefore once an index enters $Q$ via $V_{QP}$ it cannot return to $P$ until the final $V_{PQ}$ — so in the $PP$ block the **cross-coupling $V_{PQ}/V_{QP}$ appears exactly twice**, while all internal $V$'s are $Q$–$Q$ blocks and resum into

$$
\mathcal G^Q=G^R+G^RV_{QQ}G^R+G^RV_{QQ}G^RV_{QQ}G^R+\cdots
=\big[1-G^RV_{QQ}\big]^{-1}G^R,
$$

which equals $[(G^R)^{-1}-V_{QQ}]^{-1}=[\omega-H_0^Q-V_{QQ}]^{-1}$ since $(G^R)^{-1}=\omega-H_0^Q$ on $Q$.

**Physical reading.** $V_{PQ}\,\mathcal G^Q\,V_{QP}$ is the energy-dependent self-energy / optical potential that the rest bands induce on the active space. $\tilde V$ is exactly the downfolded effective potential of Feshbach projection.

> **Key point — what "second order" means.**
> The $PP$ block is *intrinsically* second order in the active–rest coupling, and **that is exact, not a truncation**. The "higher-order" coupling terms (multiple excursions in and out of $P$) are not lost: they are carried by the *active* layer, Eq. (6.1), through the resummation $[1-\tilde V G^A]^{-1}$. Thus "drop higher orders" must be stated precisely — see §5.

---

## 4. Full scattered $T$-matrix (exact reorganization)

Return to the exact equation $T=V+V(G^A+G^R)T$ and move the rest term to the left:

$$
(1-VG^R)\,T=V+VG^A T
\;\Longrightarrow\;
T=\underbrace{(1-VG^R)^{-1}V}_{T^R}+(1-VG^R)^{-1}V\,G^A T,
$$

i.e. a **renormalized Lippmann–Schwinger equation**

$$
\boxed{\;T(\omega)=T^R(\omega)+T^R(\omega)\,G^A(\omega)\,T(\omega)\;}
\qquad\Longrightarrow\qquad
T=\big[1-T^RG^A\big]^{-1}T^R=T^R\big[1-G^AT^R\big]^{-1}.
\tag{4.1}
$$

The bare potential $V$ is replaced by the rest-dressed $T^R$, and the only propagator left is the active $G^A$.

### 4.1 No double counting

Eq. (4.1) is **exact** for the original $T$. Expanding $T=T^R+T^RG^AT^R+\cdots$ and substituting $T^R=V+VG^RV+\cdots$ produces precisely the propagator strings
$$
R^{m_0}\,A\,R^{m_1}\,A\,R^{m_2}\cdots,\qquad m_i\ge 0,
$$
i.e. every sequence with arbitrarily many $G^R$ (including zero) between consecutive $G^A$. The "zero $G^R$" case is supplied by $T^R\!\to\!V$ at leading order, so adjacent $G^A$'s are also covered. This reconstructs $T=V+V(G^A+G^R)V+\cdots$ with **no omission and no repetition**. Conclusion: the $A/R$ split is just a reorganization of one and the same series; by itself it introduces no approximation.

---

## 5. Approximations enter here: rest self-dressing and the static limit

Two independent, physically transparent approximations turn the exact (4.1)/(3.2) into something cheap.

### 5.1 Coupling-second-order vs. exact ($V_{QQ}$ dressing)

Eq. (3.2) is exact. Two truncation choices exist:

- **Exact (recommended):** keep $\mathcal G^Q=[\omega-H_0^Q-V_{QQ}]^{-1}$.
- **Drop rest self-dressing ($V^{(2)}$):** set $\mathcal G^Q\to G^R$,
$$
\tilde V^{(2)}(\omega)=V_{PP}+V_{PQ}\,G^R\,V_{QP}.
\tag{5.1}
$$

The relative error of $\tilde V^{(2)}$ is $\mathcal O(\|V_{QQ}\|/\Delta)$. It is justified only when $\|V\|/\Delta\ll1$. As §7–§9 show, **keeping $V_{QQ}$ costs almost nothing** (it only changes the Sternheimer operator from $H_0^Q$ to $H_0^Q+V_{QQ}$, or equivalently adds outer iterations), so the exact form is preferred unless a deliberate leading-order estimate is wanted.

### 5.2 Static limit for the rest

Because rest bands lie far from $\omega_0$, the denominators $\omega-\varepsilon_n\;(n\in R)$ are large and nearly $\omega$-independent across the active window. Freeze them at $\omega_0$:

$$
G^R(\omega)\to G^R(\omega_0)\equiv G^R_{\rm stat},
\qquad
T^R\to \tilde V\equiv\big[1-VG^R_{\rm stat}\big]^{-1}V\quad(\omega\text{-independent}).
$$

Then **all** frequency dependence of the full $T$-matrix comes from $G^A$:

$$
\boxed{\;T(\omega)=\big[1-\tilde V\,G^A(\omega)\big]^{-1}\tilde V\;},
\qquad
T_{PP}(\omega)=\big[1-\tilde V\,G^A(\omega)\big]^{-1}\tilde V,
\tag{5.2}
$$

with $\tilde V$ evaluated once at $\omega_0$.

**Validity condition (independent of §5.1).** The static limit requires the active bandwidth to be small compared with the gap, $W_A\ll\Delta$, and $\omega_0$ taken at the window center. If $W_A$ is not small, keep the first frequency correction (§10). Note this is a *separate* condition from $\|V_{QQ}\|/\Delta\ll1$: §5.1 controls coupling strength, §5.2 controls frequency variation.

---

## 6. Assembling the full $T$-matrix (active layer)

Only the active block of (4.1) is needed, because $G^A=PG^AP$ pins both sides into $P$. Projecting (4.1) with $P$ and using $PT^R G^A = PT^RP\,G^A$ (valid since $G^A=PG^AP$):

$$
T_{PP}=\tilde V+\tilde V\,G^A\,T_{PP}
\quad\Longrightarrow\quad
\boxed{\;T_{PP}(\omega)=\big[1-\tilde V\,G^A(\omega)\big]^{-1}\tilde V\;}.
\tag{6.1}
$$

This is a **single small matrix inversion** in the active space ($\dim = N_A$, the number of active states), performed at each $\omega$ of interest, and it resums the active dynamical multiple scattering **exactly**. Whether $\tilde V$ is taken at coupling-second-order (5.1) or fully dressed (3.2), this layer is the same.

> **Do not conflate $\tilde V$ with $T_{PP}$.** $\tilde V$ — even fully converged in rest dressing — is only the *input*. The full $T$-matrix is (6.1). The two-layer picture: the **ladder of §9 feeds the all-order rest dressing into $\tilde V$ ($k$-decoupled, iterable to exact); the active dynamics are then resummed exactly by the small inversion (6.1).**

---

## 7. Sternheimer formulation of the rest dressing

To avoid explicitly enumerating the (many) high-energy rest bands in $\mathcal G^Q$, define for each active state $|n\rangle\;(n\in A)$ the **partially perturbed wavefunction**

$$
|\delta\psi_n^R(\omega_0)\rangle\equiv \mathcal G^Q(\omega_0)\,V_{QP}\,|n\rangle,
\qquad V_{QP}|n\rangle=Q\,V\,|n\rangle .
$$

It satisfies the **dressed Sternheimer equation**

$$
\boxed{\;Q\big(\omega_0-H_0-V\big)Q\;|\delta\psi_n^R\rangle=Q\,V\,|n\rangle,
\qquad P|\delta\psi_n^R\rangle=0\;}
\tag{7.1}
$$

(the **bare** version drops $V$ inside the operator: $Q(\omega_0-H_0)Q$). The effective potential matrix elements are then a small set of inner products ($m,n\in A$):

$$
\tilde V_{mn}=\langle m|V|n\rangle+\langle m|V|\delta\psi_n^R\rangle.
\tag{7.2}
$$

### 7.1 Why Sternheimer is *especially* well-behaved here

The Sternheimer operator's spectrum is $\{\omega_0-\varepsilon_k:k\in R\}$ (dressed: the $V_{QQ}$-shifted rest spectrum). Since $\omega_0$ sits inside the active window and **far** from the rest bands, these eigenvalues are all bounded away from zero: the operator is **non-singular, well-conditioned, free of small denominators**. All the dangerous near-resonant states were projected out by $Q$ and handed to the *small* active inversion (6.1). Krylov solvers (CG/GMRES) converge in a few iterations. This is the entire payoff of the $A/R$ split.

### 7.2 Cost

Only $N_A$ linear solves are needed, and — thanks to the static limit — they are solved **once** at $\omega_0$ and reused for all $\omega$. Because $V$ is a fixed difference potential, there is **no DFPT-style self-consistency cycle**: the problem is simpler than standard DFPT.

### 7.3 Hermiticity-preserving bilinear form

Instead of the one-sided inner product (7.2), use the bilinear form

$$
\langle m|V|\delta\psi_n^R\rangle
=\big\langle\delta\psi_m^R\big|\,(\omega_0-H_0^Q-V_{QQ})\,\big|\delta\psi_n^R\big\rangle
\tag{7.3}
$$

(bare version: replace $V_{QQ}\to 0$). Both sides equal the second term of $\tilde V_{mn}=(V_{PQ}\mathcal G^Q V_{QP})_{mn}$ at convergence; the right-hand side is **manifestly Hermitian** in $(m,n)$ for real $\omega_0$, so $\tilde V$ stays Hermitian even when the linear solves are not fully converged. (Derivation: $\langle\delta\psi_m^R|(\mathcal G^Q)^{-1}|\delta\psi_n^R\rangle=\langle m|V_{PQ}\mathcal G^Q V_{QP}|n\rangle$ since $(\mathcal G^Q)^{-1}|\delta\psi_n^R\rangle=V_{QP}|n\rangle$ and $\langle\delta\psi_m^R|=\langle m|V_{PQ}\mathcal G^Q$.)

---

## 8. The $k$-grid problem: why bare decouples and dressed does not

On a Bloch basis $|nk\rangle$, the decisive question is whether $V$ sits **inside the operator** or only **in the source**.

### 8.1 Bare Sternheimer is block-diagonal in $k$

$H_0$ and $\omega_0$ are lattice-periodic, hence block-diagonal in crystal momentum:

$$
Q(\omega_0-H_0)Q=\bigoplus_{k}\;Q(k)\,[\omega_0-H_0(k)]\,Q(k),
$$

and the rest projector $Q(k)=\sum_{n\in R}|nk\rangle\langle nk|$ acts $k$ by $k$ (it only selects which bands are "rest" at each $k$; it does not connect $k\to k'$). The bare equation therefore splits into $N_k$ **independent** per-$k$ solves:

$$
Q(k')\,[\omega_0-H_0(k')]\,Q(k')\;|\chi_n^{(0)}(k')\rangle=\big[V|n\rangle\big](k').
\tag{8.1}
$$

The defect potential $V$ acts **once**, only when building the source. (That single action does scatter $k\to k'$, but it is one FFT / local multiplication, $\mathcal O(N\log N)$, independent of $N_k$, and is *not* inside the iteration.) The result is the coupling-second-order $\tilde V^{(2)}$ of (5.1).

### 8.2 Dressed Sternheimer couples all $k$

A single defect breaks translational symmetry, so $V(\mathbf r)$ is **not** lattice-periodic; its Fourier components connect arbitrary $k\to k'$ and $\langle nk|V|n'k'\rangle$ is **dense** in $(k,k')$. Putting $V_{QQ}$ inside the operator (7.1) means it acts at *every* iteration, coupling all $k$. The object to be inverted is then a dense $(N_kN_b^R)\times(N_kN_b^R)$ coupled system — the $N_k\times N_k$ scaling problem. (It can also become ill-conditioned if a defect level is pushed toward $\omega_0$; §10.)

This is exactly the DFPT structure: DFPT inverts $(H_0+\alpha P-\varepsilon)$ per $(k,q)$ rather than placing $V_{\rm SCF}$ inside the operator, which is precisely why it decouples in $k$ and $q$.

---

## 9. $k$-decoupled ladder to the full rest dressing

The false dichotomy is "cheap second order vs. expensive dressed coupled solve." There is a third path: an iterative ladder that **stays $k$-decoupled at every step yet converges to the exact dressed $\tilde V$**.

### 9.1 The ladder

Neumann-expand $\mathcal G^Q$ in $V_{QQ}$; each term is built from the $k$-block-diagonal $G^R$ and a **single** action of $V_{QQ}$:

$$
|\chi_n^{(0)}\rangle=G^R\,V_{QP}|n\rangle,
\qquad
|\chi_n^{(m)}\rangle=G^R\,V_{QQ}\,|\chi_n^{(m-1)}\rangle,
\qquad
\tilde V_{mn}=\langle m|V|n\rangle+\sum_{m\ge0}\langle m|V|\chi_n^{(m)}\rangle.
\tag{9.1}
$$

- $m=0$: the coupling-second-order term — one source-$V$ action + one per-$k$ solve (8.1).
- $m=1$: third order — act with $V_{QQ}$ once on $|\chi^{(0)}\rangle$, solve bare Sternheimer again (still $k$-block-diagonal).
- and so on.

Each step inverts the **same** $k$-block-diagonal operator $A_0\equiv Q(\omega_0-H_0)Q$ (so $G^R=A_0^{-1}$ on $Q$); $V_{QQ}$ appears only in building the right-hand side (one FFT, independent of $N_k$).

### 9.2 Equivalence to the dressed Sternheimer solution

The ladder is the Richardson/Neumann iteration for (7.1). Splitting the dressed operator,

$$
\underbrace{Q(\omega_0-H_0-V)Q}_{\text{dressed}}=\underbrace{Q(\omega_0-H_0)Q}_{A_0\;(k\text{-block diagonal})}-\underbrace{V_{QQ}}_{k\text{-dense}},
$$

the partial sums $S_M=\sum_{m=0}^{M}|\chi_n^{(m)}\rangle$ satisfy

$$
A_0\,S_{M+1}=Q\,V|n\rangle+V_{QQ}\,S_M,
$$

the fixed point of which is $(A_0-V_{QQ})\,|\delta\psi_n^R\rangle=QV|n\rangle$, i.e. **exactly** the dressed Sternheimer equation. Hence

$$
S_M\;\xrightarrow{M\to\infty}\;\mathcal G^Q V_{QP}|n\rangle=|\delta\psi_n^R\rangle,
\qquad
\text{equivalently}\quad
|\chi_n\rangle=G^R\big(V_{QP}|n\rangle+V_{QQ}|\chi_n\rangle\big).
\tag{9.2}
$$

**So: yes — the order-by-order Sternheimer ladder converges to the full (all-order in $V_{QQ}$) rest dressing**, i.e. to $\tilde V=PT^RP$, *without ever leaving the $k$-decoupled per-$k$ solves.* The coupling-second-order $\tilde V^{(2)}$ is simply the $M=0$ truncation.

### 9.3 The $N_k\times N_k$ system is never formed

The dense scaling appears **only** if one explicitly builds and directly inverts the dressed operator ($\mathcal O((N_kN_b^R)^3)$, prohibitive). But the dressed system never needs an explicit matrix: both the Neumann ladder (9.1) and a Krylov method (GMRES, using $A_0$ as preconditioner and $V_{QQ}$ as an FFT matvec) require **only matrix–vector products**. Neumann is the simplest such scheme (slower); GMRES is its accelerated cousin (few steps). The real boundary is therefore **"form dense matrix and invert" (infeasible) vs. "iterative matvec" (feasible)** — and second order is just the zeroth matvec step.

---

## 10. Convergence and failure modes

**Convergence.** The ladder (and the equivalent Richardson iteration) converges geometrically with ratio

$$
\rho\big(G^R V_{QQ}\big)\sim \frac{\|V_{QQ}\|}{\Delta}.
$$

This is exactly the condition under which the $A/R$ partition is meaningful, so **whenever the separation is physically justified, the ladder converges**, often in one or two steps.

**Diagnostic (for the Anvil tests).** Monitor the successive-ratio
$$
\frac{\|\tilde V^{(m+1)}-\tilde V^{(m)}\|}{\|\tilde V^{(m)}-\tilde V^{(m-1)}\|},
$$
and the condition number / smallest eigenvalue of the per-$k$ bare operator $A_0(k)$ as a function of $V$ strength. If the ratio is not $\ll1$, or if the condition number blows up as $\|V\|$ grows, a rest state is being pushed toward $\omega_0$.

**Failure mode — deep/bound levels.** If the defect pulls a bound/resonant state out of the rest continuum toward $\omega_0$, then $\|V_{QQ}\|/\Delta\gtrsim1$ and the ladder diverges; correspondingly the dressed operator $\omega_0-H_0^Q-V_{QQ}$ approaches singularity at $\omega_0$. This is intrinsically non-perturbative. Note the *structure* (3.2) is still exact (coupling-second-order $PP$ block); it is the $V_{QQ}$ resummation that fails. Two remedies:

1. **Recommended:** move the near-resonant rest state **into the active space** (shrink the partition). The newly enlarged active block then carries that state non-perturbatively via the small dressed active solve, and the rest dressing is again weakly coupled. Physically this is usually the correct choice.
2. Directly solve the coupled dressed system (Krylov with strong preconditioning), accepting the cost.

---

## 11. Beyond the static limit: first frequency correction

If the active window is not narrow ($W_A$ not $\ll\Delta$), retain the leading frequency term:

$$
\mathcal G^Q(\omega)\approx\mathcal G^Q(\omega_0)+(\omega-\omega_0)\,\partial_\omega\mathcal G^Q\big|_{\omega_0},
\qquad
\partial_\omega\mathcal G^Q=-(\mathcal G^Q)^2.
$$

The derivative is obtained from **one additional Sternheimer solve** acting on $|\delta\psi_n^R\rangle$ (since $\partial_\omega|\delta\psi_n^R\rangle=-\mathcal G^Q|\delta\psi_n^R\rangle$), i.e.
$$
(\omega_0-H_0^Q-V_{QQ})\,|\delta\psi_n^{R,(1)}\rangle=-|\delta\psi_n^R\rangle,
$$
so the linear frequency dependence of $\tilde V$ is captured without re-solving across all $\omega$.

---

## 12. Metallic / continuum case

If the rest bands extend down toward $\omega_0$, $G^R(\omega_0)$ acquires an imaginary part and $\tilde V$ becomes **non-Hermitian** — it gains a decay width, which is the (physical) absorptive part of the optical potential. Use a complex shift $\omega_0+i\eta$: the Sternheimer operator remains non-singular but requires complex arithmetic, and the downstream active inversion (6.1) is handled with a non-Hermitian $\tilde V$. If the rest bands are genuinely gapped away from $\omega_0$, this subsection does not apply.

---

## 13. The two-layer structure at a glance

$$
\underbrace{\;V\;\xrightarrow[\text{$k$-decoupled Sternheimer ladder, §7–§9}]{\text{rest dressing, all orders in }V_{QQ}}\;\tilde V(\omega_0)\;}_{\textbf{Layer 1: rest, statically renormalized, no dense }N_k\times N_k}
\;\xrightarrow[\text{exact active resummation}]{\quad\text{Eq. (6.1)}\quad}\;
\underbrace{\;T_{PP}(\omega)=[1-\tilde V G^A(\omega)]^{-1}\tilde V\;}_{\textbf{Layer 2: active dynamics, small inversion}}
$$

- **Layer 1** is exact to second order in the active–rest coupling (Feshbach), iterable to all orders in the rest self-energy $V_{QQ}$ while remaining block-diagonal in $k$. Approximations: optional drop of $V_{QQ}$ ($\mathcal O(\|V_{QQ}\|/\Delta)$) and static limit ($\mathcal O(W_A/\Delta)$).
- **Layer 2** is an exact, $\omega$-dependent inversion in the small active space.
- Neither layer constructs or inverts a dense $(N_kN_b)\times(N_kN_b)$ matrix.

---

## 14. Implementation notes (EDI / QE / Wannier)

These notes name the concrete EDI routines and QE facilities that each step maps onto. A tag marks reuse status: **[have]** = exists in `edi-dev/src` and is directly reusable; **[QE]** = available in the linked QE libraries but not yet wired into EDI; **[new]** = must be written. The full status/dictionary/reuse audit is §16.

1. **Source $V_{QP}|n\rangle=QV|n\rangle$ — built from the parts EDI already forms.** The difference potential $V=\Delta V$ is EDI's `V_colin` (collinear) or `V_nc(:,:)` (noncollinear, scalar+B$_{xc}$), constructed in `edi.f90`/`ed_coarse.f90` by `build_vcolin_aligned` / `build_vcolin_corealign` from the supercell potentials `V_d,V_p` (`edic_mod::V_file`). **[have]** Acting $V$ on a state is exactly what `ed_coarse_full_q` already does to build $M_B$: inverse-FFT the periodic part to real space (`invfft`, `dffts%nl`, `igk_k`), multiply by `V_colin` with the supercell→primitive fold + Bloch phase, plus the Kleinman–Bylander nonlocal piece via `get_betavkb`+`calbec`+`dvan`/`dvan_so`. The rest projector $Q=1-P$ is "project out the $N_A$ Wannier-active states," the standard orthogonalization already used to disentangle the active manifold. The Sternheimer source is therefore *one reuse of the existing $V$-action followed by a projection* — no new physics. **[new]** is only the wiring into a linear solver.

2. **Per-$k$ solves (8.1).** Each $k$ is an independent solve of $Q(k)[\omega_0-H_0(k)]Q(k)$ on the **primitive** cell — embarrassingly parallel over the coarse $k$-grid, which is exactly the pool layout EDI already runs over (`mp_pools`, `poolcollect`, the panel-broadcast loop in `ed_coarse_full_q`). The operator $\omega_0-H_0(k)$ is the standard QE Hamiltonian action `h_psi` on the primitive cell already loaded by `read_file_new`; the well-conditioned, projected linear solve is a Krylov method from `KS_Solvers/libks_solvers.a` (already linked in `src/makefile`). **[QE]** No enumeration of high-energy rest bands; the single source action is $\mathcal O(N\log N)$ and independent of $N_k$.

3. **Rest dressing without a dense matrix.** Implement (9.1)/(9.2) as Neumann or GMRES with $A_0$ as preconditioner and $V_{QQ}$ as an FFT matvec (same `V_colin`-action machinery as point 1, restricted by $Q$). Truncate at the order set by the §10 diagnostic; $M=0$ recovers the coupling-second-order estimate (this is the term whose Born-level cousin EDI computes today). **[new]**

4. **Hermiticity.** Use the bilinear form (7.3), not the one-sided inner product, so $\tilde V$ stays Hermitian under incomplete linear-solver convergence (real $\omega_0$). This matches EDI's existing habit of Hermitianizing interpolated operators (`hamwan2bloch_edi`, `hamwan2bloch_with_evec` symmetrize $H(k)$ before `ZHEEV`).

5. **Frequency.** With the static limit, solve once at $\omega_0$ and reuse for all $\omega$ in (6.1); $\omega_0$ is naturally the center of EDI's transport window (`transport_win_min/max`). If the window is wide, add the single extra solve of §11 for the linear-in-$\omega$ term.

6. **Active inversion (6.1).** Small dense inversion of $[1-\tilde V G^A(\omega)]$ at each target $\omega$ (dimension $N_A=$ `nbndsub`, typically 5–11). $G^A(\omega)=\sum_{a\in A}|a\rangle\langle a|/(\omega-\varepsilon_a)$ is the explicit sum over the few active Wannier bands, whose energies $\varepsilon_a$ and eigenvectors $U(k)$ EDI already produces by interpolating `chw` (=$H_W(R)$) with `hamwan2bloch_with_evec`. The downfolded $\tilde V$ replaces today's bare $M_W(k,k')$ in the same Bloch-rotation slot ($M_B=U^\dagger\,\tilde V\,U$) used in `transport.f90` / `edmatwan2bloch_2d`, after which the rate is the usual $\propto|T_{PP}|^2\,\delta(\varepsilon)$ in place of $|M|^2$. **[new]** for $\tilde V$ and the inversion; **[have]** for the rotate-and-square-into-golden-rule tail.

---

### Pseudocode skeleton  *(annotated with EDI routine names; "→ new" = not yet in src)*

```text
# ---- Layer 1: build the downfolded effective potential V_tilde at omega_0 ----
# inputs already available in EDI: V_colin (build_vcolin_*), primitive H0 (read_file_new),
#   active manifold = nbndsub Wannier states, U(k) from hamwan2bloch_with_evec
for each active state |n>:                       # N_A = nbndsub states
    source_n = Q * V * |n>                        # reuse ed_coarse_full_q V-action:
                                                  #   invfft + V_colin fold/phase  (local)
                                                  #   get_betavkb+calbec+dvan(_so) (nonlocal)
                                                  #   then project out active manifold (Q)
    chi = solve_sternheimer(A0=Q(omega0 - H0)Q,   # → new: KS_Solvers Krylov on primitive h_psi
                            rhs=source_n,          #   per-k (mp_pools), well-conditioned
                            projector=Q)           # = |chi^(0)>  (coupling-2nd-order)

    while not converged(ratio_test):              # → new: Neumann / GMRES outer loop
        rhs   = V_QQ * chi                         # one FFT (V_colin action), no k-coupling in SOLVE
        dchi  = solve_sternheimer(A0, rhs, Q)      # same per-k operator A0
        chi  += dchi
    delta_psi[n] = chi                             # = |delta_psi_n^R>  (dressed)

for m,n in active x active:                        # Hermitian bilinear form (7.3)  → new
    V_tilde[m,n] = <m|V|n> + <delta_psi[m]| (omega0 - H0^Q - V_QQ) |delta_psi[n]>

# ---- Layer 2: exact active multiple scattering at each frequency ----
for omega in target_grid:
    GA   = sum_{a in active} |a><a| / (omega - eps_a)   # eps_a from hamwan2bloch (chw)  [have]
    Tpp  = inv( I - V_tilde @ GA ) @ V_tilde            # small nbndsub x nbndsub inversion  → new
    # then: M_B = U(ki)^dagger · Tpp · U(kf), rate ∝ |M_B|^2 δ(eps)  — transport.f90 tail  [have]
```

---

## 15. Relation to cRPA (dictionary)

| cRPA (screened interaction) | This work (multiple scattering) |
|---|---|
| bare interaction $v$ | bare defect potential $V$ |
| target/low-energy polarization $P^{A}$ | active propagator $G^A$ (kept dynamical) |
| rest polarization $P^{R}$ | rest propagator $G^R$ (kept static) |
| partially screened $W^{r}=(1-vP^R)^{-1}v$ | partial $T$-matrix $T^R=(1-VG^R)^{-1}V$ |
| fully screened $W=(1-W^rP^A)^{-1}W^r$ | full $T_{PP}=(1-\tilde V G^A)^{-1}\tilde V$ |

The structural identity is exact: rest (high-energy) channels statically renormalize the bare object; active (low-energy) channels retain the full dynamics.

---

## 16. Grounding in the EDI codebase: status, dictionary, and reuse map

This section reconciles the abstract construction above with the production EDI code at `/anvil/projects/x-che190065/rjguo/qe-7.5/edi-dev` (QE-7.5 plugin, `edi.x`). It states precisely what EDI does today, why the $A/R$ partition is *already physically present* in EDI, what maps onto what, and which pieces are reusable versus new.

### 16.1 What EDI computes today: the first Born term only

EDI is a Wannier-interpolated, supercell **difference-potential** code. Its entire output is built from the single matrix element

$$
M_{nk,mk'}=\langle\psi_{nk}\,|\,\Delta V\,|\,\psi_{mk'}\rangle,
\qquad \Delta V=V_{\rm defect}-V_{\rm pristine},
$$

with local and Kleinman–Bylander nonlocal parts, and it inserts $|M|^2$ into Fermi's golden rule. The transport kernel (`transport.f90`, `compute_transport`) is literally

```fortran
inv_tau_serta(ibnd,iki) += twopi * n_d_cell * wqf * ABS(edmatf_b(ibnd,jbnd))**2 * w_delta
```

i.e. $1/\tau_{nk}=2\pi\,n_d\,\frac1{N_k}\sum_{mk'}|M_{nk,mk'}|^2\,\delta(\varepsilon_{nk}-\varepsilon_{mk'})$ (MRTA adds the $1-\cos\theta$ factor). In the language of §1–§4 this is the **leading Born term**: $T=V+VG_0T$ truncated at $T\approx V$. EDI never applies a host resolvent $G_0$, never resums multiple scattering, and contains no $T$-matrix, Sternheimer, Feshbach, or $A/R$-projection code (verified by grep over `src/`). **Everything in §3–§13 is therefore a beyond-Born extension of the existing EDI pipeline, not a description of it.**

The current pipeline (orchestrated by `edi.f90`):

1. **Potential** — load $V_d,V_p$ (`load_supercell_pot` on-the-fly via `read_file_new`, or `load_pot_from_file` cube files), align (`build_vcolin_aligned`/`_corealign`, `pot_align='vacuum'|'core'|'none'`) → $\Delta V=$ `V_colin` (or `V_nc`).
2. **Range separation** (optional, `range_sep.f90::compute_range_separation`) — split $\Delta V=\Delta V_{\rm LR}+\Delta V_{\rm SR}$ via $\Delta\rho\!\to\!$ Coulomb kernel (2D-truncated Ismail-Beigi or 3D $8\pi/G^2$, optional Gaussian $\alpha$).
3. **Coarse $M$** (`ed_coarse.f90::ed_coarse_full_q`) — local part by `invfft`+`V_colin` fold/phase; nonlocal by `get_betavkb`+`calbec`+`dvan`/`dvan_so`; over the kept bands `nbnd_kept` (= non-`excluded_band`).
4. **Bloch→Wannier** (`edbloch2wan.f90::edbloch2wane`) — double FT to $M(R_e,R_p)$; **diagonal/low-rank** structure exploited and SVD-compressed (`tt_compress.f90::svd_compress_edmat`, <5% singular values → <1% error).
5. **Wannier→Bloch** (`wan2bloch.f90`) — $H_W=$`chw` interpolated to $\varepsilon(k),U(k)$ (`hamwan2bloch_with_evec`); $M$ interpolated (`edmatwan2bloch_2d`) and rotated $M_B=U^\dagger(k_i)M_W U(k_f)$; long-range $M^{\rm LR}$ added semi-analytically (`lr_matelem.f90::compute_mlr_at_kpair`, EPW dipole/Berry-connection form).
6. **Transport** (`transport.f90`, `delta_weights.f90`, `bz_symmetry.f90`) — golden-rule rates, IBZ symmetry, SERTA/MRTA mobility, Fermi-level bisection.

### 16.2 The key bridge: EDI's Wannier window *is* the active space $A$

The single most important observation for an implementer: **EDI already partitions the band structure into exactly the $A/R$ split this note needs, for an unrelated (Wannierization) reason.** Disentanglement keeps `nbndsub` Wannier functions spanning the bands near $E_F$ inside the windows `dis_win_min/max` (and frozen `dis_froz_min/max`), excluding deep/core bands (`bands_skipped`→`exclude_bands`) and high-energy bands. Therefore:

- **Active $A$ ($P$)** $\equiv$ the `nbndsub` kept Wannier bands (`n_wannier`, `num_bands`, `u_mat`, `wann_centers`; kept set = non-`excluded_band`, `nbndep`/`ibndkept`). These are the strongly-coupled, dynamical, near-$E_F$ states — precisely where transport lives (`transport_win_min/max`).
- **Rest $R$ ($Q$)** $\equiv$ the bands EDI throws away from the Wannier manifold: the deep states removed by `exclude_bands` and the high-energy states above the outer window. These are the "distant bands, treated statically" of §2.
- **$G^A(\omega)$** is built from data EDI already has: the interpolated active eigenvalues $\varepsilon_a(k)$ from `chw`/`hamwan2bloch`, as $\sum_{a}|a\rangle\langle a|/(\omega-\varepsilon_a)$.
- The energy gap $\Delta$ of §10 is the disentanglement-window-to-rest-band separation; the static-limit validity $W_A\ll\Delta$ (§5.2) is the statement that the Wannier window is narrow compared with its distance to the excluded bands.

So the $A/R$ machinery is conceptually free in EDI — it is the Wannier window with a sign flip ($Q=1-P$). What is *missing* is using $R$ dynamically (the Sternheimer dressing) and resumming $A$ (the active inversion).

### 16.3 Dictionary: theory object → EDI object → source

| This note | EDI / QE object | Where |
|---|---|---|
| Defect potential $V=\Delta V$ | `V_colin` / `V_nc` (real-space $\Delta V$ on SC grid) | `edic_mod.f90`; built in `edi.f90`, `ed_coarse::build_vcolin_*` |
| Action of $V$ on a state | `invfft`+`V_colin` fold/phase (local) + `get_betavkb`/`calbec`/`dvan(_so)` (nonlocal) | `ed_coarse.f90::ed_coarse_full_q`, `get_betavkb.f90` |
| First Born element $\langle\psi|V|\psi'\rangle=T^{(1)}$ | `edmatkq`/`edmatf_b` ($M_B(k,k')$) | `ed_coarse_full_q`, `transport.f90` |
| Host $H_0$, $G_0=(\omega-H_0)^{-1}$ | primitive-cell QE Hamiltonian (`read_file_new`, `h_psi`); not assembled as a resolvent | QE `PW/src` |
| Active projector $P$, space $A$ | kept Wannier manifold `nbndsub` (non-`excluded_band`) | `wann_common.f90`, `global_var.f90` |
| Rest projector $Q$, space $R$ | excluded + high-energy bands ($Q=1-P$) | (implicit; **new**) |
| Active propagator $G^A(\omega)$ | $\sum_a|a\rangle\langle a|/(\omega-\varepsilon_a)$ from $\varepsilon_a,U(k)$ | `wan2bloch.f90::hamwan2bloch_with_evec`, `chw` |
| Eigvecs $U(k)$ for the Bloch rotation | `evec`/`U(k)` (filukk `u_kc`) | `read_filukk_edi`, `hamwan2bloch_with_evec` |
| Downfolded effective potential $\tilde V$ | would replace bare $M_W$ in the rotate-and-square slot | (**new**; slot is `transport.f90`/`edmatwan2bloch_2d`) |
| Partially perturbed $|\delta\psi_n^R\rangle$ (7.1) | rest-projected Sternheimer solution | (**new**; solver from `KS_Solvers`) |
| Full active block $T_{PP}(\omega)$ (6.1) | $[1-\tilde V G^A]^{-1}\tilde V$, dim `nbndsub` | (**new**) |
| Golden-rule rate $\propto|T_{PP}|^2\delta$ | `inv_tau_serta/mrta` loop | `transport.f90` |

### 16.4 Reuse map — what exists vs. what must be built

**Reusable as-is [have]:**
- Construction and alignment of $\Delta V$ (`build_vcolin_*`, `range_sep`), and the *action* of $\Delta V$ on Bloch states (local FFT fold + KB nonlocal) inside `ed_coarse_full_q` — this is the Sternheimer **source** generator.
- Active eigensystem $\varepsilon_a(k),U(k)$ from `chw` via `hamwan2bloch_with_evec`; the Bloch rotation $M_B=U^\dagger(k_i)\,\cdot\,U(k_f)$; the golden-rule/SERTA/MRTA tail (`transport.f90`). $\tilde V$ or $T_{PP}$ drops straight into the slot now occupied by $M_W$/$M_B$.
- Pool parallelism over $k$ (`mp_pools`, `poolcollect`, panel broadcast) — the natural home for the $k$-decoupled per-$k$ solves of §8/§9.

**Available but not yet wired [QE]:**
- Krylov linear solvers (`KS_Solvers/libks_solvers.a`, already in `src/makefile`) and the primitive-cell `h_psi` operator — the engine for the per-$k$ Sternheimer solve (8.1). Note EDI does **not** link `LR_Modules`, so QE's `solve_linter` is not directly present; the §7 solver must be assembled from `h_psi`+`KS_Solvers` (it is simpler than DFPT — see §7.2, no self-consistency since $\Delta V$ is fixed).

**Must be written [new]:**
1. The rest projector $Q=1-P$ as an "orthogonalize against the `nbndsub` active states" operation, and the source $QV|n\rangle$.
2. The per-$k$ Sternheimer solve (7.1)/(8.1) with $Q(\omega_0-H_0)Q$ (bare) and the §9 Neumann/GMRES ladder for the $V_{QQ}$ dressing.
3. Assembly of $\tilde V_{mn}$ via the Hermitian bilinear form (7.3).
4. The small active inversion (6.1) per target $\omega$, replacing $M_W$ by $T_{PP}$ in the existing rotate-and-square path.

### 16.5 Long-range / metallic caveats meet EDI's range separation and low-rank structure

- **§12 (metallic/continuum, non-Hermitian $\tilde V$, complex $\omega_0+i\eta$) and the long-range tail.** EDI already faces the long-range problem head-on: a charged or polar defect makes $\Delta V\sim1/r$, which is *not* short-ranged, breaks the locality that makes $M(R,R')$ compact, and is split off as $\Delta V_{\rm LR}$ (`range_sep.f90`) and treated semi-analytically (`lr_matelem.f90`, EPW dipole form). For the $T$-matrix construction the same split is the natural remedy: apply the Sternheimer/$A/R$ machinery to the **short-range** $\Delta V_{\rm SR}$ (where $\|V\|/\Delta\ll1$ and the ladder of §9 converges), and keep the long-range Coulomb channel in the existing analytic LR path. If active bands extend to $\omega_0$ (metallic host), §12's complex shift $\omega_0+i\eta$ applies and $\tilde V$ acquires the physical absorptive width.
- **Low-rank / SVD.** EDI finds $M(R,R')$ empirically low-rank (`tt_compress.f90`: <5% of singular values for <1% error; see `plan_zSVD_update.md`). This is the multiple-scattering echo of why the *active layer* is cheap: the physically relevant scattering lives in a small subspace (the `nbndsub` active block), exactly the space in which $\tilde V$ and the (6.1) inversion are formed. The SVD low rank is empirical evidence that $N_A$ is small and the §6 inversion is genuinely tiny.

### 16.6 Minimal implementation path inside EDI

A staged route that reuses the most and adds the least:

1. **Coupling-second-order $\tilde V^{(2)}$, bare ($V_{QQ}\!=\!0$), static.** For each active Wannier state build $QV|n\rangle$ with the existing `ed_coarse_full_q` $V$-action, solve (8.1) per $k$ with a `KS_Solvers` Krylov method on the primitive `h_psi`, assemble $\tilde V^{(2)}_{mn}$ via (7.3). This is the smallest delta over today's code (the source machinery is shared) and already goes beyond Born for the active dynamics once fed through step 3.
2. **Dress the rest** with the §9 ladder (`V_QQ` matvec via the same `V_colin` action) until the §10 ratio test passes — full all-orders $\tilde V$, still $k$-decoupled.
3. **Active resummation:** at each transport-window $\omega$ compute $T_{PP}=[1-\tilde V G^A(\omega)]^{-1}\tilde V$ (dim `nbndsub`), rotate $M_B=U^\dagger\,T_{PP}\,U$, and feed the existing `inv_tau` loop. Validate against the Born limit by checking $T_{PP}\to\tilde V\to V_{PP}$ when the dressing and the active resummation are switched off (must reproduce today's $M$).
4. **Frequency/range-sep refinements:** add the §11 single extra solve if $W_A\not\ll\Delta$; route the long-range channel through `range_sep`/`lr_matelem` (16.5).

The validation hooks of §10 (successive-ratio test; smallest eigenvalue / condition number of $A_0(k)$ vs. $\|V\|$) should be emitted on `ionode` in the same diagnostic style EDI already uses (e.g. the Wannier-interpolation error check and the range-separation roundtrip check), so a deep/bound-level failure mode is caught early.

---

*End of note.*
