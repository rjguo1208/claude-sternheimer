# The rest-space Sternheimer ladder: why second order is not enough, and solving the higher-order dressing term by term

EDT's active/rest downfolding writes the dressed potential as $\tilde V(\omega)=M+\Sigma(\omega)$: $M$ is the bare coupling within the active subspace, and $\Sigma$ folds the virtual transitions through the rest subspace (the high-energy bands that were folded out) back into the active space. The $\Sigma$ currently in the code goes only to **second order (Born)** — it uses the bare host propagator $G^0_{QQ}$ and discards the multiple scattering of the defect potential inside the rest space, $\Delta V_{QQ}$. On the MoS$_2$ S-vacancy this approximation **over-screens** the conduction-derived $e$ level by about $1.1$ eV (see Figs 13–15 on the [results page](results.html#sec-3)).

This page does three things: (i) derives the exact Feshbach rest self-energy; (ii) explains **why second order is not a good approximation**; and (iii) gives the complete formulae and implementation plan for the **order-by-order Sternheimer ladder** (including how the nonlocal part attaches). Notation: projectors $P$ (active) and $Q=1-P$ (rest); host Hamiltonian $H_0$, defect potential $\Delta V=V_{\rm defect}-V_{\rm host}$; block notation $\Delta V_{PQ}\equiv P\,\Delta V\,Q$ and so on.

## 1. Notation and the exact downfolding (Feshbach)

The host is diagonal in its Bloch eigenbasis, $H_0|n\mathbf k\rangle=\varepsilon_{n\mathbf k}|n\mathbf k\rangle$. Active and rest are **both** spanned by host Bloch states, so $H_0$ does **not** couple the two blocks:
$$P H_0 Q=0\quad(\text{Bloch states are eigenstates of }H_0),$$
and the $P$–$Q$ coupling comes **only** from $\Delta V$. Write the active-block bare coupling as $M\equiv\Delta V_{PP}$.

Solving the stationary equation $(H_0+\Delta V-\omega)|\psi\rangle=0$, substituting $|\psi\rangle=|\psi_P\rangle+|\psi_Q\rangle$ and eliminating $|\psi_Q\rangle$ from the $Q$-component equation:
$$|\psi_Q\rangle=(\omega-H_{QQ})^{-1}\,\Delta V_{QP}\,|\psi_P\rangle,\qquad
H_{QQ}\equiv Q(H_0+\Delta V)Q=H_0^{QQ}+\Delta V_{QQ},$$
gives the effective Hamiltonian in the active subspace
$$\boxed{\,H_{\rm eff}(\omega)=P H_0 P+\underbrace{M+\Sigma(\omega)}_{\tilde V(\omega)}\,},\qquad
\Sigma(\omega)=\Delta V_{PQ}\,\big(\omega-H_0^{QQ}-\Delta V_{QQ}\big)^{-1}\,\Delta V_{QP}. \tag{1.1}$$

$\Sigma(\omega)$ is the rest dressing. **Note that the denominator contains $\Delta V_{QQ}$** — the multiple scattering of the defect potential within the rest subspace, which the exact $\Sigma$ keeps to **all orders**. This $\tilde V$ is precisely the input to the active-space $T$-matrix $T_{PP}=[1-\tilde V G^A]^{-1}\tilde V$, so the accuracy of $\Sigma$ determines the whole pipeline. The defect levels are the self-consistent eigenvalues of $\det[\omega-H_{\rm eff}(\omega)]=0$.

## 2. Why the second-order (Born) approximation over-screens

The current code drops $\Delta V_{QQ}$ from the denominator and uses the bare host propagator
$$G^0_{QQ}(\omega)\equiv\big(\omega-H_0^{QQ}\big)^{-1}=Q(\omega-H_0)^{-1}Q,\qquad
\Sigma^{(2)}(\omega)=\Delta V_{PQ}\,G^0_{QQ}(\omega)\,\Delta V_{QP}. \tag{2.1}$$

Dyson-expanding the exact resolvent of (1.1),
$$\big(\omega-H_0^{QQ}-\Delta V_{QQ}\big)^{-1}=\sum_{p\ge0}\big(G^0_{QQ}\,\Delta V_{QQ}\big)^p G^0_{QQ},$$
gives
$$\Sigma(\omega)=\sum_{p\ge0}\Sigma^{(2+p)},\qquad
\Sigma^{(2+p)}=\Delta V_{PQ}\,G^0_{QQ}\big(\Delta V_{QQ}\,G^0_{QQ}\big)^p\,\Delta V_{QP}. \tag{2.2}$$
Second order keeps only $p=0$, throwing away every rest-internal multiple scattering that contains $\Delta V_{QQ}$ ($p\ge1$).

**(a) Convergence.** The series (2.2) converges $\iff$ the spectral radius
$$\rho\equiv\rho\!\big(G^0_{QQ}(\omega)\,\Delta V_{QQ}\big)<1,\qquad
\frac{\lVert\Sigma^{(m+1)}\rVert}{\lVert\Sigma^{(m)}\rVert}\xrightarrow{\,m\to\infty\,}\rho.$$
Second order is only trustworthy when $\rho\ll1$. The S-vacancy is a **deep potential** (the defect cube well is $\sim-28$ Ry), so $\Delta V_{QQ}$ is strong and $\rho\sim\mathcal O(1)$ — the series is far from converged, and truncating at second order can be off by an $\mathcal O(1)$ factor.

**(b) Physical picture: the intermediate propagation is unscreened.** Using the resolvent identity $G_{QQ}=G^0_{QQ}+G^0_{QQ}\Delta V_{QQ}G_{QQ}$, second order's error is
$$\Sigma-\Sigma^{(2)}=\Delta V_{PQ}\,G^0_{QQ}\,\Delta V_{QQ}\,G_{QQ}\,\Delta V_{QP},$$
which is **first order** in $\Delta V_{QQ}$ and not small. $\Sigma^{(2)}$ treats the intermediate rest states as propagating like **free host bands** ($G^0$); but those rest states are themselves strongly scattered by the deep defect potential $\Delta V_{QQ}$. For the conduction-derived $e$ (the antibonding combination of the three Mo dangling bonds), its second-order coupling to the rest conduction manifold ($\varepsilon_q>\varepsilon_e$),
$$\Sigma^{(2)}_{ee}=\sum_{q\in Q}\frac{|\langle e|\Delta V|q\rangle|^2}{\varepsilon_e-\varepsilon_q}<0$$
is strongly attractive (negative denominators). With a bare $G^0$ that attraction is **unscreened** and pulls $e$ far too deep; putting $\Delta V_{QQ}$ back means letting the intermediate rest states rearrange in the defect field, and the higher-order terms partly cancel that bare attraction.

**(c) Comparison with the data.** Explicit puts these conduction bands directly into $P$ with no rest dressing at all, and gives $e$ at $+1.35$ eV above the VBM (Fig 14); second-order dressed gives $+0.36$ eV (Fig 15). In other words, **the all-order resummation has to lift $e$ from $+0.36$ back to $\sim+1.2\text{–}1.35$** (DFT gives $+1.19$) — the higher-order terms net out about $1$ eV of excess binding. That is the quantitative meaning of "second order over-screens".

> A second and **independent** approximation stacks on top of this: the code evaluates $\Sigma$ at a single static reference $\omega_0=\varepsilon_{\rm VBM}$ rather than self-consistently at each level's own energy. For $e$, which is far from $\omega_0$, even the second-order term is evaluated at the wrong frequency. The ladder on this page addresses the **truncation order**; the frequency-dependence layer requires either evaluating $\Sigma(\omega)$ at the level energy or sweeping the full $\Sigma_{\rm rest}(\omega)$.

## 3. The Sternheimer ladder: order-by-order derivation

Sternheimer's central idea is to **avoid the explicit rest-band sum** and solve a (projected) linear equation instead. For each active source $|b\rangle$, define the first-order rest response
$$|\chi^0_b\rangle\equiv G^0_{QQ}\,\Delta V_{QP}\,|b\rangle
\;\Longleftrightarrow\;
Q(\omega-H_0)Q\,|\chi^0_b\rangle=Q\,\Delta V\,|b\rangle,\quad|\chi^0_b\rangle\in Q. \tag{3.1}$$
(In QE this is solved by projected CG, with $H_0$ supplied by `h_psi`.) The second-order self-energy is then $\Sigma^{(2)}_{ab}=\langle a|\Delta V|\chi^0_b\rangle$.

**Recursion for the higher-order responses.** Peeling $\big(G^0_{QQ}\Delta V_{QQ}\big)^p$ in (2.2) apart layer by layer, define
$$|\chi^p_b\rangle\equiv\big(G^0_{QQ}\,\Delta V_{QQ}\big)^p|\chi^0_b\rangle
\;\Longleftrightarrow\;
Q(\omega-H_0)Q\,|\chi^p_b\rangle=Q\,\Delta V\,Q\,|\chi^{p-1}_b\rangle,\quad p\ge1. \tag{3.2}$$
**The key point**: every order uses **the same** operator $Q(\omega-H_0)Q$ and only changes the right-hand side (scattering the previous order's response once more by $\Delta V_{QQ}$). The CG setup is reused verbatim, and each extra order costs only "build one right-hand side + one solve".

**Expressing each order of the self-energy through the responses (symmetric pairing).** On the real axis $G^0_{QQ}$ is Hermitian, so $\langle a|\Delta V_{PQ}G^0_{QQ}=\langle\chi^0_a|$; splitting the chain of (2.2) at any intermediate point,
$$\boxed{\;\Sigma^{(2)}_{ab}=\langle a|\Delta V|\chi^0_b\rangle,\qquad
\Sigma^{(2+p)}_{ab}=\langle\chi^i_a|\,\Delta V_{QQ}\,|\chi^j_b\rangle,\quad i+j=p-1\ \ (p\ge1).\;} \tag{3.3}$$
(With a broadening $\eta$, $G^0$ is non-Hermitian; use the asymmetric form $\Sigma^{(2+p)}_{ab}=\langle a|\Delta V|\chi^p_b\rangle$, or solve the left and right responses separately.)

**Cost structure: one solve buys two orders, and the odd orders come free.** By (3.3), having the responses $\chi^0,\ldots,\chi^{S-1}$ ($S$ solves) gives exactness to order $2S+1$:

| responses solved | number of solves $S$ | self-energy order reachable | the "free" term from that solve |
|---|---|---|---|
| $\chi^0$ | $1$ | $\Sigma^{(2)},\ \Sigma^{(3)}$ | $\Sigma^{(3)}=\langle\chi^0|\Delta V_{QQ}|\chi^0\rangle$ |
| $+\,\chi^1$ | $2$ | $\Sigma^{(4)},\ \Sigma^{(5)}$ | $\Sigma^{(5)}=\langle\chi^1|\Delta V_{QQ}|\chi^1\rangle$ |
| $+\,\chi^2$ | $3$ | $\Sigma^{(6)},\ \Sigma^{(7)}$ | $\Sigma^{(7)}=\langle\chi^2|\Delta V_{QQ}|\chi^2\rangle$ |
| $\ \ \vdots$ | $S$ | exact to $\Sigma^{(2S+1)}$ | $\langle\chi^{S-1}|\Delta V_{QQ}|\chi^{S-1}\rangle$ |

In particular, **third order $\Sigma^{(3)}=\langle\chi^0_a|\Delta V_{QQ}|\chi^0_b\rangle$ needs no new solve at all** — it is one extra contraction of the second-order responses already computed, i.e. free. That gives the cheapest possible reliability test for second order: if $r_3\equiv\lVert\Sigma^{(3)}\rVert/\lVert\Sigma^{(2)}\rVert$ is already not small, second order cannot be trusted.

**The all-order endpoint (resummation).** $S\to\infty$ is equivalent to solving the complete Feshbach equation with $\Delta V_{QQ}$ included directly,
$$Q(\omega-H_0-\Delta V)Q\,|X_b\rangle=Q\,\Delta V\,|b\rangle,\qquad \Sigma_{ab}=\langle a|\Delta V|X_b\rangle, \tag{3.4}$$
i.e. putting $\Delta V_{QQ}$ into the CG operator and solving to all orders in one go (full derivation, matvec and complex solver in §6). The value of the order-by-order ladder is that it **expands this non-perturbative solution into a monitorable sequence**: it can be truncated on demand when $\rho<1$, and the magnitudes of $\Sigma^{(3)},\Sigma^{(5)},\ldots$ read off $\rho$ directly and tell you whether second order deserves trust.

## 4. How the nonlocal part attaches, and the implementation plan

The right-hand side of the recursion (3.2) is $Q\,\Delta V\,Q\,|\chi^{p-1}\rangle$ — the defect potential applied to a rest wavefunction. $\Delta V$ splits into local and nonlocal parts:
$$\Delta V=\Delta V_{\rm loc}(\mathbf r)+\Delta V_{\rm nl},\qquad
\Delta V_{\rm nl}=\sum_{IJ}|\beta_I\rangle\,\Delta D_{IJ}\,\langle\beta_J|.$$
$\Delta V_{\rm nl}$ is KB **separable** (for a vacancy it is essentially minus the removed atom's projector); $\Delta V_{\rm loc}$ is the aligned supercell local potential difference (the cube). Applying $\Delta V$ to $|\chi^{p-1}\rangle$:

```text
given the rest response |chi^{p-1}>:
1. local:     FFT |chi^{p-1}> to the supercell real-space grid -> multiply pointwise by dV_loc(r) -> FFT back to G space
2. nonlocal:  for each projector J compute c_J = <beta_J | chi^{p-1}>
              assemble sum_{I,J} |beta_I> dD_{IJ} c_J            (separable, cheap)
3. add -> dV|chi^{p-1}>
4. project onto Q:  |R> = (1-P) dV|chi^{p-1}>
                        = dV|chi^{p-1}> - sum_{a in P} |a><a| dV|chi^{p-1}>
5. solve Sternheimer:  Q(w - H0)Q |chi^p> = |R>     (projected CG, reusing second order's operator)
```

Two implementation essentials:

- **The local part couples the $\mathbf k$ channels.** The second-order solve (3.1) is **decoupled** in $\mathbf k$ (at fixed source, the right-hand side is independent per $\mathbf k$). But the defect breaks translational symmetry, and $\Delta V_{\rm loc}$ connects different $\mathbf k$ through the supercell FFT: the higher-order right-hand side $\Delta V Q|\chi^{p-1}\rangle$ scatters one $\mathbf k$'s response into others. So **third order and above are no longer per-$\mathbf k$ independent solves**, and this supercell FFT field machinery is needed — the main new engineering effort. The nonlocal part remains a separable projector sum throughout: cheap and $\mathbf k$-decoupled.
- **Reuse the right-hand sides for the contractions.** In (3.3), $\Sigma^{(2+p)}_{ab}=\langle\chi^i_a|\Delta V_{QQ}|\chi^j_b\rangle=\langle\chi^i_a|\big(\Delta V|\chi^j_b\rangle\big)$, and $\Delta V|\chi^j_b\rangle$ is exactly the right-hand-side vector already built when computing $\chi^{j+1}$. So the self-energy contractions cost no extra $\Delta V$ application, just vector inner products.

## 5. Convergence criteria, the correspondence with explicit calculations, and next steps

**How to use the ladder to judge second order.** The practical order of operations: compute $\Sigma^{(3)}$ for free first, then look at the ratio $r_3=\lVert\Sigma^{(3)}\rVert/\lVert\Sigma^{(2)}\rVert$.

| $r_3$ | meaning | what to do |
|---|---|---|
| $\ll1$ | second order is trustworthy | $\tilde V\approx M+\Sigma^{(2)}$, done |
| $\lesssim1$ | slow convergence | climb order by order ($+1$ solve to fifth order, then look at $r_5$) |
| $\gtrsim1$ | Born **diverges** | going order by order cannot save it — switch to the all-order Feshbach (3.4), or move these bands into the active space and treat them explicitly |

**Measured $r_3$ (using the explicit $M$, exact, with the full $k$ coupling).** Taking $P$ = bands 7–17 (the block's active window), $Q$ = bands 18–28 (the near-gap conduction rest), and the physical potential $W=M/N_k$, direct matrix algebra gives $\Sigma^{(2)}=W_{PQ}\,g\,W_{QP}$ and $\Sigma^{(3)}=W_{PQ}\,g\,W_{QQ}\,g\,W_{QP}$:

| $\omega$ | $\lVert\Sigma^{(2)}\rVert_F$ | $\lVert\Sigma^{(3)}\rVert_F$ | $r_3$ | $\rho(W_{QQ}g)$ |
|---|---|---|---|---|
| $\omega_0=\varepsilon_{\rm VBM}$ (the block's static reference) | $8.72$ eV | $3.22$ eV | $\mathbf{0.37}$ | $0.60$ |
| $\omega\approx\varepsilon_e$ ($e$'s own energy) | $10.35$ eV | $4.64$ eV | $\mathbf{0.45}$ | $0.70$ |

$\rho<1$ but $r_3\approx0.4$ — the Born series **converges, but very slowly** (third order still carries $40\%$ of second order), landing squarely in the **"$\lesssim1$, keep climbing"** row above: second order is not trustworthy.

**The ladder really does converge, and it can be verified self-consistently.** On the same data, the $e$ level (relative to the VBM, static $\Sigma(\omega_0)$) as a function of order:

| treatment | $e$ (relative to VBM) |
|---|---|
| bare $M$ (P=7–17, no rest) | $+1.49$ |
| $+\,\Sigma^{(2)}(\omega_0)$ | $+1.30$ (second-order overshoot) |
| $+\,\Sigma^{(2)}+\Sigma^{(3)}(\omega_0)$ | $+1.40$ (the free third order pulls it back) |
| $+\,\Sigma_{\rm full}(\omega_0)$ (all-order resummation) | $+1.37$ |
| self-consistent Feshbach (all-order, self-consistent $\omega$) | $\mathbf{+1.35}$ |
| explicit all-band 7–28 (target) | $+1.35$ ✓ |

The self-consistent Feshbach reproduces explicit's $+1.35$ **exactly**, verifying that the downfolding framework is self-consistent; and the free third-order term is measured to pull the second-order overshoot from $+1.30$ back to $+1.40$, exactly the order-by-order convergence one expects.

**A note on the $Q$ truncation.** Here $Q$ extends only to band 28 (the ceiling of the explicit data), which is the dominant near-gap rest, and second order overshoots it only **mildly** ($+1.30$ vs $+1.35$). But the block uses the **full rest, 18–150**: the extra high bands push $\rho$ toward $1$, which is what produced the drastic $+0.36$ over-screening seen earlier ($e$ dropping $\sim1.1$ eV from the bare $+1.49$). So **the measured $r_3\approx0.4$ and $\rho\approx0.6$–$0.7$ are a lower bound on the full problem** — the full rest can only be closer to divergence and needier of a ladder or of the all-order Feshbach. (Measurement script: `edt/run/ladder_r3.py`.)

**Measured in the code (in-code cross-validation).** The free $\Sigma^{(3)}$ term has landed in the EDT block code (the `do_sigma3` switch, single-state mode): retain the per-channel $\chi^0$ (rather than discarding it after use), then do one cross-channel double sum over $\Delta V$, $\Sigma^{(3)}_{aa}=\sum_{k',k''}\langle\chi^0(k')|\Delta V(k'\!\leftarrow\!k'')|\chi^0(k'')\rangle$ — **with no new CG solve**. For band 14, $k=1$ (full rest 18–150):

| quantity | Fortran (full rest) | explicit ($Q=18$–$28$) |
|---|---|---|
| $\Sigma^{(2)}_{aa}$ self-consistency check | reconstruction $=$ the code's stored `Sgblk`, identical to 7 digits ✓ | — |
| $\Sigma^{(2)}_{aa}$ | $-0.045$ eV | $-0.014$ eV |
| $\Sigma^{(3)}_{aa}$ | $+0.070$ eV | $+0.0053$ eV |
| $r_{3,aa}$ | $\mathbf{1.54}$ | $0.37$ |

**(1)** The $\Sigma^{(2)}_{aa}$ self-consistency check **passes exactly** (reconstructing from the saved $\chi^0$ $=$ the code's stored `Sgblk`), proving the cross-channel machinery is correct. **(2)** The full rest gives $r_{3,aa}=1.54>1$ — **third order is larger than second, so the Born series diverges**, landing in the **"$\gtrsim1$, divergent"** row of the criterion table above. This **upgrades "the full rest is closer to divergence" from an inference to a measurement**: the truncated rest ($r_3=0.37$) is a lower bound, and the full rest really does cross $1$ (the distant conduction bands make $\Sigma^{(3)}$, with its two rest sums, grow much faster than $\Sigma^{(2)}$: measured, $|\Sigma^{(3)}|$ is 13× larger while $|\Sigma^{(2)}|$ is only 3× larger). So the final conclusion is firm: **the full rest cannot be saved order by order and requires either the all-order Feshbach (resolvent inversion, not a series) or explicit treatment**; what diverges is the **order-by-order** expansion, and the explicit self-consistent Feshbach still gives exactly $+1.35$. (Caveat: the intermediate $\Delta V_{QQ}$ takes only the local part, and this is the single-state diagonal $r_{3,aa}$; but $1.54\gg1$ is robust. Switch `do_sigma3`, job `edt/run/edt_sigma3.{in,slurm}`.)

**Two routes to the same physics.** The same $e$ level can be computed two ways: **(i)** put the conduction bands into $Q$ and take the dressing to all orders (this page's ladder / the all-order Feshbach); **(ii)** put them into $P$ and diagonalize explicitly with no dressing (the 21-band explicit on the results page). Given enough bands the two converge to the same answer. Explicit already gives a clean $a_1\oplus e$ ($e$ at $+1.35$, Fig 14), so **the ladder's convergence target is known** and can be used to validate the implementation.

**All-order closure (measured, direct resolvent).** Using the explicit-60 matrices, the all-order dressing for $Q=$ bands 18–70 ($N_Q{=}7174$, a band-converged rest) was computed by **direct inversion** (one eigendecomposition of $H_{QQ}$, no iteration): the $e$ level (relative to the VBM) order by order = bare $+1.495$ → $\Sigma^{(2)}$ **$+0.732$** (a disaster, with $a_1$ squeezed out of the gap) → $\Sigma^{(2)}{+}\Sigma^{(3)}$ **$+1.595$** (the divergent series **oscillating**, overcorrecting) → all-order static $+1.223$ → **self-consistent-$\omega$ all-order $+1.205$ = the explicit-60 all-band value, agreeing to the digit** (DFT $+1.19$). The order-by-order ratios for this rest are $r_3{=}0.74,\ r_4{=}0.91,\ r_5{=}1.21$ — the ratios cross $1$, completing the divergence scale (between $0.37$@18–28 and $1.54$@full rest). The conclusion is reinforced: **order by order is unsalvageable, and all-order resummation lands it in one step**; details and spectral functions are on the [results page §3](results.html#sec-3) (Fig 17). Also: the measured spectrum of $H_{QQ}$ comes no closer than $5.5$ eV to $\omega_0$ (a well-conditioned static resolvent) — so the earlier MINRES stagnation was **not** caused by a near-null mode, but by the mismatch between the kinetic-energy Jacobi preconditioner and the $+11.8$ Ry core region (so when reviving the iterative line, change the preconditioner first).

**Implementation status.** The free $\Sigma^{(3)}$ term is in place (`do_sigma3`). The first two MVP steps — ① reuse the second-order solve to obtain $\chi^0$; ② compute $\Sigma^{(3)}$ for free in-code and measure $r_3$ — are **complete** (both the explicit route's $r_3=0.37$–$0.45$ and the in-code full-rest $r_{3,aa}=1.54$ have been measured). Because $r_3>1$ (the full rest diverges), the third step falls straight through to **"switch to the all-order Feshbach (3.4) / explicit"**: the order-by-order ladder does not converge for this system, so the production route is all-order resolvent inversion or moving the bands into the active space and treating them explicitly. The acceptance criterion (lifting the dressed $e$ from second order's $+0.36$ eV back to explicit/DFT's $\sim+1.2\text{–}1.35$ eV) has been met by the 21-band explicit (Fig 14).

## 6. All-order Feshbach Sternheimer: the production route once order-by-order diverges

§5 measured the order-by-order ladder to diverge on the full rest ($r_3=1.54$), so second order's over-screening **cannot** be repaired by adding a few orders. The way out is the all-order Feshbach that **does not expand but inverts directly** — still a Sternheimer linear solve, only with $\Delta V_{QQ}$ inside the operator. This is the complete derivation of eq. (3.4).

**The linear equation.** Substituting $H_{QQ}=Q(H_0+\Delta V)Q=QHQ$ back into (1.1):
$$\Sigma(\omega)=\Delta V_{PQ}\,(\omega-QHQ)^{-1}\,\Delta V_{QP}.$$
For each active source $|b\rangle$ define the all-order rest response $|X_b\rangle\equiv(\omega-QHQ)^{-1}\Delta V_{QP}|b\rangle\in Q$; multiplying on the left by $(\omega-QHQ)$ and using $Q|X_b\rangle=|X_b\rangle$ (so $QHQ|X_b\rangle=QH|X_b\rangle$):
$$\boxed{\,Q(\omega-H_0-\Delta V)Q\,|X_b\rangle=Q\,\Delta V\,|b\rangle\,},\qquad
\Sigma_{ab}=\langle a|\Delta V|X_b\rangle=\langle S_a|X_b\rangle, \tag{6.1}$$
where $|S_a\rangle=Q\Delta V|a\rangle$ is **the same** source as in second order; symmetrically, $\Sigma=S^\dagger A^{-1}S$ with $A=Q(\omega-H_0-\Delta V)Q$. **Exact, all-order, untruncated.**

**The only difference from the second-order Sternheimer: one extra $\Delta V$ in the matvec.**

| | second order (Born) | all-order Feshbach |
|---|---|---|
| operator $A$ | $Q(\omega-H_0)Q$ | $Q(\omega-H_0-\Delta V)Q$ |
| matvec | `h_psi` $+\,Q$ | `h_psi` $-\,\Delta V+Q$ |
| $\mathbf k$ structure | block diagonal → decoupled per $\mathbf k$ | $\Delta V$ couples them → one large equation over all $\mathbf k$ |

```text
A|v> ,  v in Q   (each CG iteration):
  1. t = omega*v - h_psi(v)      # (omega - H0) v       [per-k, already exists]
  2. t = t - dV(v)               # subtract dV v : local supercell FFT (build_V_folded, cross-channel) + separable nonlocal
  3. t = Q t  (apply_Qproj)      # project back into rest
```

That is, §4's $\Delta V$ application routine (the same cross-channel machinery as `do_sigma3`) goes from "compute $\Sigma^{(3)}$ once" to "**called once per matvec**".

**Why this converges when the order-by-order expansion does not — inversion $\neq$ series.** Write $A=A_0-B$ with $A_0=Q(\omega-H_0)Q$ and $B=Q\Delta V Q$. Order by order is the Neumann series $A^{-1}=\sum_{p\ge0}A_0^{-1}(BA_0^{-1})^p$, which needs $\rho(A_0^{-1}B)=\rho(G^0\Delta V_{QQ})<1$ — measured at $\sim1.54$, divergent. **But $A^{-1}$ itself exists**; the convergence of Krylov/CG applied directly to $AX=S$ depends only on the spectrum of $A$, not on that $\rho$. By analogy: at $x=1.54$ the series $1+x+x^2+\cdots$ diverges while $\tfrac{1}{1-x}=-1.85$ is perfectly well defined — **what diverges is the expansion, not the number**. With $\eta>0$, $A=Q(\omega+i\eta-H)Q$ is always invertible.

**Conditioning and solvers.** With the rest lying above $\omega_0$, $A_0\succ0$; after adding $-\Delta V_{QQ}$, $A$ remains positive definite provided no rest state is pulled below $\omega$ (any defect bound state should already be inside $P$) → real CG. If a rest resonance approaches $\omega$ ($A$ indefinite or ill-conditioned) → go to $\omega+i\eta$ with a complex solver (the code already has `rest_split='complex'` / `ccgsolve_all`), with $\eta$ doubling as a regularizer; the active space is pinned into $Q$ with $\alpha P$ deflation.

**$\omega$ self-consistency → exact levels.** $\Sigma(\omega)$ depends on $\omega$; the defect levels are the self-consistent roots of $\det[\omega-H_0^{PP}-M-\Sigma(\omega)]=0$, obtained by iterating $\omega$. At self-consistent $\omega$, the all-order Feshbach gives **exact** levels — the same quantity that reproduced explicit's $+1.35$ in §5.

**Cost and where it sits.** One **all-$\mathbf k$** Sternheimer solve per source, $n_{\rm iter}$ steps with one $\Delta V$ application per step (the cost scale of $\Sigma^{(3)}$), over all $N_A$ sources of the block — and it **converges** (the series does not). It and the 21-band explicit corroborate each other: one puts conduction in $Q$ and dresses to all orders, the other puts it in $P$ and does not dress — and since order-by-order diverges, these two are the routes to take.

**Cost estimate (anchored to measurements, 36 ranks / 1 node).**

| run | size | measured WALL |
|---|---|---|
| second-order block (`h_psi` matvec, per $\mathbf k$) | $N_A=1584$, full BZ | $1$ h $59$ m |
| explicit 21-band $M$ (born_only, no solve) | $N_A=3002$ | $46$ m |
| single-state $\Sigma^{(3)}$ (`do_sigma3`, folded $\Delta V$) | $1$ source $+$ one cross-channel pass | $10$ m |

The all-order Feshbach $\approx$ second order's `h_psi` part ($\sim2$ h) $+\ N_A\times n_{\rm iter}\times t_{\Delta V}$. All of the cost sits in that per-matvec cross-channel $\Delta V$, and $t_{\Delta V}$ **depends entirely on the implementation**:

| $\Delta V$ implementation | one $t_{\Delta V}$ (single source, all channels) | after $N_A\,n_{\rm iter}\!\sim\!10^5$ |
|---|---|---|
| folded double sum (`do_sigma3` as it stands) | $\sim5$ min ($144^2$ supercell foldings) | **months, infeasible** |
| supercell-FFT (needs writing; one supercell FFT covers all channels) | $\sim1$–$4$ s | $\Delta V$ part $\sim1$–$2$ h, **feasible** |

A gap of **$10^2$–$10^3\times$**: the folded route (the one used for the single-state $\Sigma^{(3)}$) recomputes $144^2$ supercell foldings at every inner step, which inside an iteration means months. The correct approach applies $\Delta V$ with **one supercell FFT** (assemble $X$'s $144$ channel components into supercell real space $\to$ multiply by $\Delta V(\mathbf r)$ $\to$ reassemble, covering all $\mathbf k$ in one pass). Via supercell-FFT: **$\sim3$–$4$ h per $\omega$** (indefinite/ill-conditioned operators requiring BiCGStab plus more iterations add another $\times2$–$3$), and **self-consistent $\omega$ (the defect levels, $\sim5$–$10$ iterations) $\sim1$–$3$ days per node**.

**For comparison: explicit is $1$–$2$ orders of magnitude cheaper.** Explicit 21-band (measured at **$46$ min, one pass**) gives the same $e=+1.35$ with no iteration and no $\Delta V$-in-matvec — $\sim10$–$50\times$ cheaper even than an optimized all-order Feshbach. So for **defect levels**, explicit is the practical production route; the all-order Feshbach is irreplaceable only when (i) the full frequency-dependent $\Sigma(\omega)$ is needed (a spectral-function self-energy), or (ii) the active space cannot be enlarged any further and the bands must stay in $Q$ and be dressed to all orders — and even then it must go via supercell-FFT, with `do_sigma3`'s folded route fit only for **one-off** diagnostics. The concrete code implementation plan (static-$\omega_0$, a cube-anchored `apply_dV`, an element-by-element acceptance ladder, a source-parallel layout, and the post-mortem of the v1 acceptance) is in the [Feshbach implementation plan (v2)](feshbach-implementation.html).
