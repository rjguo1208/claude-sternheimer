# Deflated ladder: explicit all-order rest ($\le$150) + Sternheimer tail ($>$150)

## 1. Motivation: why neither pure scheme survives

Stage A of the Sdisp probe (one S displaced $+0.30$ Å in the 6$\times$6 MoS$_2$ supercell,
non-SOC, truth = supercell diagonalization) left two hard facts on the table:

1. **The bare (2nd-order) Sternheimer rest over-binds.** With the 11-band active manifold and
   the rest folded at $\omega_0$ through the *pristine* propagator, a spurious pair of
   conduction-derived bound states appears in the gap ($-103/-38$ meV below the CBM; the true
   gap is clean), robust against every static-fold variant (both anchors, raised window,
   semicore-in-active). Dressing the $\le$150-band rest to all orders in $\Delta V$ removes the
   pair and lands the CBM doublet at $+4.6/+13.6$ meV vs the exact $+6.5$ meV — but the
   $>$150 tail, re-added at 2nd order, brings the artifact back at half depth. The tail
   (42% of the total rest dressing by the T2 gate) must be dressed too.
2. **The plain Neumann ladder for the full rest diverges.** Power iteration on the explicit
   $\le$150 block gives $\rho(G^0 V_{RR}) = 1.215$ (cross-checked on two machines), despite a
   5.46 eV spectral gap and $|V_{RR}|_{\max} = 36$ meV per element: the localized $\Delta V$
   couples the $2\times10^4$ rest states coherently — $\lVert\Sigma_{\rm tail}^{AA}\rVert_2 =
   7.5$ eV against 29 meV single elements. Same physics as the classic divergence of the Born
   series for potentials strong enough to bind.

The cure is to split the rest and let each half be handled by the method that is exact for it:
the dangerous, coherent low part by an **explicit all-order inverse**, and only the remote tail
by iteration.

## 2. Space partition and notation

Within the plane-wave Hilbert space at the wavefunction cutoff (measured $n_{\rm pw}\simeq
24{,}205$ per $k$; 144 $k$-points):

$$\mathcal{H}=A\;\oplus\;R_1\;\oplus\;R_2$$

| block | content | projector | dimension |
|---|---|---|---|
| $A$ | active manifold, bands 7–17 | $P_A$ | $11\times144=1584$ |
| $R_1$ | explicit rest: semicore 1–6 + conduction 18–150 | $P_1$ | $139\times144=20{,}016$ |
| $R_2$ | tail: everything the NSCF never represented, $P_2=1-P_{150}$ | $P_2$ | $\sim 3.4\times10^6$ (implicit) |

Write $V\equiv\Delta V$, $V_{XY}=P_X V P_Y$, $H_1=P_1H_0P_1=\mathrm{diag}(\varepsilon_i)$,
$H_2^0=P_2H_0P_2$. The object to compute is the Feshbach rest self-energy with $V$ to all
orders inside the rest:

$$\Sigma(\omega_0)=V_{AR}\,\big[\omega_0-H_R^0-V_{RR}\big]^{-1}\,V_{RA},\qquad R=R_1\oplus R_2 .$$

## 3. Schur elimination of $R_1$ (exact)

On $R_1\oplus R_2$ the resolvent to invert is the block operator

$$\omega_0-H_R=\begin{pmatrix}\omega_0-H_1-V_{11} & -V_{12}\\[2pt] -V_{21} & \omega_0-H_2^0-V_{22}\end{pmatrix}.$$

Define the two solvable resolvents:

$$D_1=(\omega_0-H_1-V_{11})^{-1}=W\,(\omega_0-\mu)^{-1}W^\dagger
\qquad\text{(explicit: one } 20{,}016^2 \text{ eigendecomposition, in hand)}$$

$$D_2=P_2\,(\omega_0-H_0)^{-1}P_2
\qquad\text{(implicit: Sternheimer CG with all 150 explicit bands lifted by the }\alpha\text{-shift)}$$

Standard Schur inversion then gives, with no approximation,

$$\boxed{\;\Sigma(\omega_0)=\underbrace{V_{A1}\,D_1\,V_{1A}}_{\Sigma^{R_1}_{\rm dressed}}
\;+\;\langle\tilde s_b|\,G_{22}\,|\tilde s_a\rangle\;}$$

where both tail-channel objects are *dressed by* $D_1$:

$$\tilde s_a=P_2\,V\big(1+D_1P_1V\big)\,|\psi_a\rangle
\qquad\text{(entrance channel: virtual } R_1 \text{ excursion before entering the tail)}$$

$$G_{22}=\big[\omega_0-H_2^0-W_{22}\big]^{-1},\qquad
W_{22}=\underbrace{V_{22}}_{\text{direct tail scattering}}
+\underbrace{V_{21}\,D_1\,V_{12}}_{R_1\text{-mediated round trip}} .$$

The first term is exactly the object already computed in the explicit-dressing verdict run.
All coherent low-energy scattering — the driver of the $\rho=1.215$ divergence — is locked
inside the exact $D_1$; only vertices that genuinely pass through the tail remain for
iteration.

## 4. The tail ladder (the only iterative part)

For $x_a\equiv G_{22}\tilde s_a$, iterate

$$x^{(0)}=D_2\,\tilde s_a,\qquad
\boxed{\;x^{(n+1)}=D_2\big[\tilde s_a+W_{22}\,x^{(n)}\big]\;},\qquad
\Sigma^{\rm tail}_{ba}=\langle\tilde s_b|x_a\rangle .$$

Every ingredient is an executable operation: $W_{22}$ = real-space fold-multiply by $\Delta V$
+ projections against the stored 150-band wavefunctions (ZGEMM) + small dense $D_1$ algebra in
the $\mu$-eigenbasis; $D_2$ = one CG solve. Note that the $D_2$ solves are *fast*: the tail
spectrum starts $\gtrsim 40$ eV above $\omega_0$, so the preconditioned CG condition number is
$\sim$17 (10–15 iterations), far better than the bare-arm solves with their 5.4 eV nearest
neighbor. The equivalent GMRES form solves $(1-D_2W_{22})\,x=D_2\tilde s$ with the same
operator applications and typically fewer iterations.

## 5. Measured convergence (the certificate)

The contraction rate was measured directly by power iteration on $G_{\rm pre}W$ using the
explicit all-band block, with a *model tail* $R_2'=$ bands $N_1{+}1..150$ (fully explicit, so
$\rho$ is exact) and the real construction mimicked block by block:

| split $N_1$ | model tail starts (eV) | $\rho$ | verdict |
|---|---|---|---|
| — (plain ladder, no deflation) | 5.46 | **1.215** | diverges |
| 90 | 25.7 | 0.574 | converges, slow |
| 120 | 32.8 | 0.323 | good |
| 140 | 37.5 | **0.153** | fast |

The decay is monotone and accelerating; the physical $>$150 tail sits farther still, giving
the extrapolated working point $\rho\approx0.1$–$0.2$: **4–5 rungs** (GMRES $\sim$3–4) for
$10^{-3}$ accuracy, i.e. a production cost of $\approx$4–5$\times$ the bare Sternheimer arm.

## 6. Built-in limits and validation gates

| truncation of the full formula | reduces to | gate |
|---|---|---|
| $D_1\to$ bare, $W_{22}\to0$, $\tilde s\to P_2V\psi$ | bare Sternheimer arm | bit-level regression |
| keep $D_1$, drop $W_{22}$ | explicit dressing + 2nd-order tail | must land between the V3/V2 brackets |
| full formula at a model split $N_1<150$ | fully explicit exact answer | T2-dressed comparison |

Physical acceptance: the spurious CBM-side pair must vanish and the CBM doublet land at the
exact $+6.5$ meV, with the $+215/+231$ meV resonances recovered simultaneously (the bare arm
could only get one side right per anchor).

## 7. Implementation notes (design frozen, not yet built)

1. The CG projector shift must lift **all 150 explicit bands** (the current code lifts only the
   11 active ones) — otherwise $D_2$ leaks back into $R_1$.
2. $\tilde s_a$ is built once per source and reused across rungs; all $D_1$ algebra is dense
   GEMM work in the 20,016-dimensional $\mu$-basis.
3. $W_{22}$ mixes all 144 $k$ (the dilute single-defect $\Delta V$ has continuous Fourier
   support — the commensurability-class block structure applies only to the supercell-array
   problem, as verified numerically: cross-class $|M|$ up to 0.51 vs 0.20 within class).
   The multi-$k$ vectors are gathered in source batches.
4. Convergence dial: raising the NSCF band count (150 $\to$ 300, $\sim$40 min) pushes the
   $R_1/R_2$ boundary up and lowers $\rho$ further.
