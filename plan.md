# plan.md — **EDT**: Electron–Defect *T*-matrix via downfolding + Sternheimer

A QE-based package that computes the electron–defect **$T$-matrix** beyond the first
Born approximation, using the active/rest (cRPA-style) downfolding and the
$k$-decoupled Sternheimer solve described in [`research.md`](research.md). It is built
as a Quantum ESPRESSO plug-in that **reuses the EDI code**
(`/anvil/projects/x-che190065/rjguo/qe-7.5/edi-dev`) for everything to do with the
difference potential, the Wannier rotation, and the matrix-element / transport
machinery; the genuinely new code is the Sternheimer rest-dressing layer and the
small active-space resummation.

> Working name: **EDT** (Electron–Defect *T*-matrix), executable `edt.x`. It can ship
> either as a sibling plug-in next to `edi-dev/` or as a new run-mode inside EDI; this
> plan assumes a sibling that links EDI's compiled objects (see §13).

---

## 0. The algorithm in one paragraph (and where each piece comes from)

For each active Bloch state $|n k_i\rangle$ we (i) form the source $|s\rangle = Q\,\Delta V\,|n k_i\rangle$
[**EDI** $\Delta V$-action], (ii) solve the per-$k'$ Sternheimer equation
$Q(\omega_0-H_0)Q\,|\chi\rangle = |s\rangle$ for the rest-space response **on a dense full-BZ $k'$-grid** (the
defect scatters into a continuum of $k'$, not the supercell folding set — §5) [**NEW**, QE `ccgsolve_all`+`h_psi`],
(iii) optionally dress it to all orders in $V_{QQ}$ by a $k$-decoupled ladder [**NEW**],
then (iv) assemble the downfolded effective potential
$\tilde V_{m k_f,\,n k_i}=\langle m k_f|\Delta V|n k_i\rangle+\langle\chi_{m k_f}|(\omega_0-H_0^Q-V_{QQ})|\chi_{n k_i}\rangle$
[bilinear form, Hermitian]. Finally (v) the active-space dynamics are resummed by the small inversion
$T_{PP}(\omega)=[1-\tilde V\,G^A(\omega)]^{-1}\tilde V$ [**NEW**, small dense solve], (vi) $\tilde V$/$T_{PP}$
are Wannierized and interpolated to a fine grid [**EDI** `edbloch2wane`/`edmatwan2bloch_2d`], and
$|T_{PP}|^2$ replaces $|M|^2$ in the golden-rule rate [**EDI** `compute_transport`].

Setting the rest response to zero recovers $\tilde V\to\Delta V$ and, with the active
resummation off, $T_{PP}\to M$ — i.e. **exactly today's EDI first-Born result**. That is the
primary regression test (§14).

---

## 1. Theory recap (see `research.md` for derivations)

Active projector $P$ (bands near $\omega_0$), rest $Q=1-P$, $[P,H_0]=0$. Host propagator splits
$G_0=G^A+G^R$. Two layers:

- **Layer 1 (rest, static):** Feshbach downfolding gives the *exact* (2nd-order-in-coupling) effective potential
  $$\tilde V(\omega)=V_{PP}+V_{PQ}\,\mathcal G^Q(\omega)\,V_{QP},\qquad \mathcal G^Q=[\omega-H_0^Q-V_{QQ}]^{-1},$$
  obtained from the **dressed Sternheimer equation**
  $$Q(\omega_0-H_0-V)Q\,|\delta\psi_n^R\rangle = Q\,V\,|n\rangle,\qquad P|\delta\psi_n^R\rangle=0 ,$$
  whose operator is lattice-periodic ⇒ **block-diagonal in $k$** (bare) and reachable to all orders in
  $V_{QQ}$ by a $k$-decoupled Neumann/GMRES ladder (`research.md` §7–§9).
- **Layer 2 (active, dynamic):** exact resummation
  $$T_{PP}(\omega)=[1-\tilde V\,G^A(\omega)]^{-1}\tilde V .$$

EDI today computes only $M=\langle\psi|\Delta V|\psi'\rangle=T^{(1)}$ (first Born). EDT adds Layers 1–2.

---

## 2. Architecture

```
                ┌─────────────────────────  EDT pipeline  ─────────────────────────┐
 QE NSCF  ─┐    │ A. inputs        B. source        C. Sternheimer    D. dress      │
 (primit.) │    │  ΔV=V_colin       |s>=QΔV|nk_i>     Q(ω0−H0)Q|χ>=|s>   ladder V_QQ  │
 SC pots  ─┼──▶ │  U(k),H_W=chw  ─▶ (EDI ΔV-action)─▶ (QE ccgsolve)  ─▶ (Neumann)─┐  │
 (def/pri) │    │  active set                                                    │  │
           │    │ E. assemble Ṽ  ◀───────────────────────────────────────────────┘  │
           │    │  Ṽ = <m|ΔV|n> + <χ_m|(ω0−H0^Q−V_QQ)|χ_n>   (Hermitian bilinear)    │
           │    │ F. active T_PP = [1−Ṽ G^A(ω)]^{-1} Ṽ   (small dense solve)         │
           │    │ G. Wannierize Ṽ/T_PP → M(R,R') → fine grid → transport |T_PP|^2   │
           │    └───────────────────────────────────────────────────────────────────┘
           └─ reuse EDI: load_supercell_pot, build_vcolin_*, range_sep, get_betavkb,
              read_filukk_edi, read_hr_file, hamwan2bloch_with_evec, edbloch2wane,
              edmatwan2bloch_2d, compute_transport
```

### Module map

| Module / file | Status | Role |
|---|---|---|
| `edt.f90` | **new** | main program: input, orchestration of stages A–G |
| `edt_input.f90` | **new** | `&edt_nml` namelist (extends `edinput_nml`) |
| `edt_partition.f90` | **new** | active/rest definition $P(k),Q(k)$; `apply_Qproj` |
| `edt_source.f90` | **new** | $|s\rangle=Q\,\Delta V|nk_i\rangle$ (wraps EDI ΔV-action as a *ket*) |
| `edt_sternheimer.f90` | **new** | `ch_psi_rest`, `solve_sternheimer_k`, complex shift |
| `edt_dress.f90` | **new** | $k$-decoupled $V_{QQ}$ ladder (Neumann/GMRES) |
| `edt_vtilde.f90` | **new** | assemble $\tilde V_{mk_f,nk_i}$ (bilinear form, Hermitian) |
| `edt_active.f90` | **new** | $G^A(\omega)$ + small inversion $T_{PP}$ |
| `ed_coarse.f90` | **reuse** | `load_supercell_pot`,`build_vcolin_*`,`read_filukk_edi`,`read_edmatw_2d_file`, the ΔV local-fold + KB nonlocal kernels |
| `range_sep.f90` | **reuse** | LR/SR split (route SR through Sternheimer, LR analytic) |
| `get_betavkb.f90` | **reuse** | supercell KB projectors $|\beta^{d/p}(k)\rangle$ |
| `wan2bloch.f90` | **reuse** | `hamwan2bloch_with_evec` → $\varepsilon_a(k),U(k)$ for $G^A$ |
| `edbloch2wan.f90` | **reuse** | $\tilde V(k_i,k_f)\to\tilde V(R,R')$ (double FT) |
| `transport.f90` | **reuse** | golden-rule rate / SERTA-MRTA with $|T_{PP}|^2$ |
| `edic_mod.f90` | **reuse** | `V_file`, `V_d/V_p`, `V_colin`, shared state |
| QE `LR_Modules` (`liblrmod.a`) | **link** | `cgsolve_all` / `ccgsolve_all`, `cg_psi` |
| QE `PW/src` (`libpw.a`) | **link** | `h_psi`, `g_psi`, `calbec`, FFT, `read_file_new` |

---

## 3. Active/rest partition in code

The decisive design choice (`research.md` §8): the projector lives **in the source/solution**, never
inside a $k$-coupled operator, so each $k$ is independent.

- **Active set** $A(k)$: kept bands ( `excluded_band` filter ) whose interpolated energy lies in the
  window `[active_win_min, active_win_max]` (= EDI's transport window / Wannier disentanglement window).
  Stored as a per-$k$ logical mask `is_active(ibnd,ik)` plus the count `n_act(ik)`.
- **$\omega_0$** = window center (input `omega0`, default = midpoint or VBM/$E_F$).
- **$Q$** = "project out the active manifold at this $k$." Implemented as Gram–Schmidt against the active
  Bloch states (norm-conserving ⇒ $S=1$), exactly the standard DFPT `orthogonalize`/`P_c^+` step.

```fortran
SUBROUTINE apply_Qproj(npw, npol, ik, nact, evc_act, psi)
  !! In-place: psi <- Q(ik) psi = psi - sum_{a in active(ik)} |a><a|psi>
  !! evc_act(npwx*npol, nact) holds the active Bloch states at ik (periodic part).
  USE kinds,    ONLY : dp
  USE mp,       ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  IMPLICIT NONE
  INTEGER,     INTENT(IN)    :: npw, npol, ik, nact
  COMPLEX(dp), INTENT(IN)    :: evc_act(:, :)
  COMPLEX(dp), INTENT(INOUT) :: psi(:)
  COMPLEX(dp) :: ovlp(nact)
  INTEGER     :: a
  ! ovlp(a) = <evc_act(a)|psi>   (sum over G, reduce over band-group)
  CALL ZGEMV('C', npw*npol, nact, (1._dp,0._dp), evc_act, SIZE(evc_act,1), &
             psi, 1, (0._dp,0._dp), ovlp, 1)
  CALL mp_sum(ovlp, intra_bgrp_comm)
  DO a = 1, nact
     psi(1:npw*npol) = psi(1:npw*npol) - ovlp(a) * evc_act(1:npw*npol, a)
  END DO
END SUBROUTINE apply_Qproj
```

---

## 4. Stage A — inputs (reuse EDI verbatim)

Nothing new. EDT calls the exact EDI entry points:

```fortran
CALL read_file()                                   ! primitive NSCF (QE) -> et, xk, evc, igk_k
CALL load_supercell_pot(defect_prefix , defect_outdir , V_d)   ! ed_coarse.f90
CALL load_supercell_pot(pristine_prefix, pristine_outdir, V_p)
CALL build_vcolin_aligned(V_d, V_p, V_colin, SIZE(V_colin), vac_shift)   ! ΔV on SC grid
IF (range_sep) CALL compute_range_separation(...)              ! ΔV = ΔV_SR + ΔV_LR
CALL read_filukk_edi(filukk, nbnd, nkstot, nks, nbndsub)       ! U(k) (u_kc)
CALL read_hr_file(prefix, nbndsub, nrr, ndegen, irvec, chw)    ! H_W(R) = chw
```

The active eigensystem $\{\varepsilon_a(k),U(k)\}$ comes from `hamwan2bloch_with_evec` (used both to
define $A(k)$ and to build $G^A$ in Stage F). **Run EDT on $\Delta V_{\rm SR}$** (the short-range part)
so $\lVert V\rVert/\Delta\ll1$ and the rest ladder converges (§7.4); keep the long-range Coulomb channel
in EDI's analytic `lr_matelem` path (`research.md` §16.5).

---

## 5. Stage B — Sternheimer source $|s\rangle = Q\,\Delta V\,|n k_i\rangle$ (reuse EDI ΔV-action)

> **The rest grid must span *all* $k$ — not the supercell-commensurate set.**
> The supercell is only a device to *capture the isolated, all-space defect potential* $\Delta V$. The
> physical defect breaks translation symmetry completely and scatters $k_i$ into a **continuum** of $k'$,
> so the rest-space resolvent
> $G^R(\omega_0)=\sum_{n\in R}\!\int_{\rm BZ}\!\frac{dk'\,|nk'\rangle\langle nk'|}{\omega_0-\varepsilon_{nk'}}$
> is a **BZ integral**: the source $|s\rangle$ and the Sternheimer response must be sampled on a **dense,
> full-BZ rest $k'$-grid**, *not* on the $N_{\rm sc}$ supercell-commensurate channels $k_i\oplus g$.
> We obtain $\Delta V$'s action at an *arbitrary* momentum transfer $q=k'-k_i$ from the **localized
> real-space** $\Delta V$ via EDI's exact-$q$ fold `build_V_folded(q)` — valid because $\Delta V$ has
> decayed to $\approx 0$ inside the supercell, the very same mechanism EDI uses for arbitrary fine-grid
> matrix elements (`ed_direct_from_files`). The rest-grid density is an explicit **convergence parameter**
> `rest_nk*` (baseline = the primitive NSCF grid, densified until $\lVert\tilde V\rVert$ converges; cf. the
> codex "Rest k-grid convergence" 4→12→144 sweep). The operator stays $k'$-block-diagonal, so this is just
> "more independent per-$k'$ solves," never a $k$-coupled system.

The only change from EDI is that we keep the **ket** $\Delta V|nk_i\rangle$ resolved per output channel
$k'$ (over the full-BZ rest grid above), instead of contracting it into a matrix element. EDI already
builds the folded potential `V_folded` for arbitrary $q$ and the KB projectors; we reuse them.

### 5.1 Local part (real-space fold of `V_colin`, then FFT to the $k'$ channel)

In `ed_coarse_full_q`, the local matrix element is
`mlocal = Σ_r conj(u_{k_i}) · V_folded^{q}(r) · u_{k_f}(r) / nnr` with
`q = k_f − k_i`. The corresponding **ket** in channel $k'$ is `ζ(r) = V_folded^{q=k'−k_i}(r) · u_{n k_i}(r)`,
then FFT-forward onto the $k'$ plane-wave set.

```fortran
SUBROUTINE dV_local_ket(ik_i, ibnd, ik_f, zeta_g)
  !! zeta_g(1:npw(ik_f)) = [ΔV_loc |ψ_{ibnd,ik_i}>] projected on channel k_f.
  USE fft_base,       ONLY : dffts
  USE fft_interfaces, ONLY : invfft, fwfft
  USE klist,          ONLY : igk_k, ngk, xk_cryst => xk   ! crystal coords used for q-phase
  USE edic_mod,       ONLY : V_colin, V_d
  IMPLICIT NONE
  INTEGER,     INTENT(IN)  :: ik_i, ibnd, ik_f
  COMPLEX(dp), INTENT(OUT) :: zeta_g(:)
  COMPLEX(dp), ALLOCATABLE :: psir(:), Vfold(:)
  ! 1) u_{ibnd,ik_i}(r): scatter evc(:,ibnd) onto FFT grid via dffts%nl(igk_k(:,ik_i)), invfft
  ! 2) Vfold(r) = build_V_folded(q = xk_cryst(:,ik_f) - xk_cryst(:,ik_i))   ! EXACT q (ed_coarse.f90 logic)
  ! 3) psir(r) <- Vfold(r) * psir(r)
  ! 4) fwfft, gather coefficients at the k_f G-set: zeta_g(ig) = psir(dffts%nl(igk_k(ig,ik_f)))
  ! (npol>1: repeat per spinor component)
END SUBROUTINE dV_local_ket
```

`build_V_folded(q)` is lifted directly from `ed_coarse_full_q` (the `V_folded_kf` fold loop with the
`d1/d2/d3` crystal-coordinate phase — the *exact* $q=k_f-k_i$, **no BZ wrapping**, per EDI's
`double_ft_bug` fix).

### 5.2 Nonlocal part (KB difference, reuse `get_betavkb`)

$\Delta V_{\rm NL}|\psi\rangle=\big(V^{d}_{\rm NL}-V^{p}_{\rm NL}\big)|\psi\rangle
=\sum_{a,ij}|\beta^{d}_i(k_f)\rangle D^{d}_{ij}\langle\beta^{d}_j(k_i)|\psi\rangle-(d\!\to\!p)$.

```fortran
SUBROUTINE dV_nonlocal_ket(ik_i, ibnd, ik_f, zeta_g)
  USE uspp,    ONLY : dvan, dvan_so, nh => nh
  USE becmod,  ONLY : calbec
  IMPLICIT NONE
  ! vkb_d_ki = get_betavkb(ngk(ik_i), igk_k(:,ik_i), xk(:,ik_i), ..., V_d%nat, V_d%tau, V_d%ityp)
  ! bec_d(jkb) = <beta^d_j(k_i)|psi_{ibnd,ik_i}>            (calbec)
  ! vkb_d_kf = get_betavkb(ngk(ik_f), igk_k(:,ik_f), xk(:,ik_f), ..., V_d%nat, ...)
  ! coeff_d(ikb) = sum_jkb dvan(ih,jh,nt) * bec_d(jkb)      ! same (na,nt,ih,jh) bookkeeping as ed_coarse
  ! zeta_g += sum_ikb vkb_d_kf(:,ikb) * coeff_d(ikb)
  ! ... subtract the pristine analog (V_p) ...
END SUBROUTINE dV_nonlocal_ket
```

### 5.3 Assemble the source over all channels and project onto rest

```fortran
DO ikp = 1, nrest                         ! channels k' over the FULL-BZ rest grid (rest_nk*), not k_i(+)g
   CALL dV_local_ket   (ik_i, ibnd, ikp, s(:,ikp))   ! build_V_folded uses the EXACT q = k'(ikp) - k_i
   CALL dV_nonlocal_ket(ik_i, ibnd, ikp, tmp); s(:,ikp) = s(:,ikp) + tmp
   CALL apply_Qproj(ngk(ikp), npol, ikp, n_act(ikp), evc_act(:,:,ikp), s(:,ikp))
END DO
```

---

## 6. Stage C — the per-$k$ Sternheimer solve (NEW core)

> **What Stage C and Stage D are — and what "second order" means here.**
> Stage C is the **bare** per-$k'$ Sternheimer solve; its output (assembled in Stage E) is the
> **coupling-second-order** effective potential $\tilde V^{(2)}=V_{PP}+V_{PQ}\,G^R\,V_{QP}$, with
> $G^R=[\omega_0-H_0^Q]^{-1}$. Two things that name does **not** mean:
>
> - **"Second order" is the *exact* order in the active–rest *coupling* $V_{PQ}$ (Feshbach), not a
>   truncation.** The $PP$ block intrinsically contains $V_{PQ}/V_{QP}$ exactly twice (`research.md` §3.1);
>   higher *coupling* orders are not lost — they are restored by Layer 2's resummation
>   $T_{PP}=[1-\tilde V G^A]^{-1}\tilde V$ (Stage F).
> - **The solve itself is exact.** The linear solver inverts the full $G^R$ (all rest bands at once via
>   `h_psi`), so Stage C is *not* a low-order / approximate solver. Its **only** approximation is dropping
>   the rest self-dressing $V_{QQ}$.
>
> **Stage D removes exactly that one approximation** — it is the $k$-decoupled ladder that adds $V_{QQ}$ to
> *all* orders, $\mathcal G^Q=[\omega_0-H_0^Q-V_{QQ}]^{-1}$. **Stage C is just the $m=0$ rung of Stage D:**
> each higher rung acts once with $V_{QQ}$ and **re-solves the same bare operator** $A_0=Q(\omega_0-H_0)Q$
> (only the right-hand side changes). So there are three *independent* "orders":
>
> | knob | Stage C (bare, $m{=}0$) | Stage D (ladder, $m{\ge}1$) | higher orders restored by |
> |---|---|---|---|
> | active–rest coupling $V_{PQ}$ | exactly 2 — *exact* (Feshbach) | exactly 2 | Layer 2 / $T_{PP}$ (Stage F) |
> | rest self-energy $V_{QQ}$ | 0 (dropped) | all orders | Stage D |
> | rest-band sum in $G^R$ | all bands — *exact* | all bands | (already exact) |

We solve, **independently for each channel $k'$ on the full-BZ rest grid** (§5 callout; `research.md` §8.1),
$$Q(k')\big(\omega_0 - H_0(k')\big)Q(k')\,|\chi(k')\rangle = |s(k')\rangle .$$

**Indefiniteness / arbitrary $\omega_0$.** On $Q$ the operator $\omega_0-H_0$ is *indefinite* (rest bands lie
both below and above $\omega_0$), so plain real CG is not appropriate. Use a **complex shift**
$\omega_0\to\omega_0+i\eta$ and QE's complex solver `ccgsolve_all` (this also covers the metallic /
continuum case, `research.md` §12). For a strictly gapped problem one may instead split $R$ into
$R^{\pm}$ (below/above $\omega_0$) and use the real `cgsolve_all` on each definite block.

### 6.1 The $A$-operator wrapper (analog of QE `ch_psi_all`)

```fortran
SUBROUTINE ch_psi_rest(n, psi, A_psi, e, ik, m)
  !! A_psi = (H0(ik) - e) psi + alpha * P_act psi      [e = omega0, REAL]
  !! ccgsolve_all carries the complex shift omega0 -> omega0 + i*eta via its freq_c arg,
  !! so ch_psi only needs the real reference e = omega0.
  !! The alpha*P_act term lifts the active eigenvalues so A is nonsingular on the
  !! full space while the iterate stays in the rest space (cf. ch_psi_all's P_v).
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER,     INTENT(IN)  :: n, ik, m
  COMPLEX(dp), INTENT(IN)  :: psi(:,:), e
  COMPLEX(dp), INTENT(OUT) :: A_psi(:,:)
  COMPLEX(dp), ALLOCATABLE :: hpsi(:,:), ppsi(:,:)
  INTEGER :: i
  CALL h_psi(npwx, n, m, psi, hpsi)                 ! QE: H0(ik) psi  (NC: no S)
  DO i = 1, m
     A_psi(1:n,i) = hpsi(1:n,i) - e*psi(1:n,i)
  END DO
  ! + alpha * P_act psi   (project onto active and add back, alpha ~ 2*(rest bandwidth))
  CALL build_Pact_psi(n, ik, m, psi, ppsi)
  A_psi(1:n,1:m) = A_psi(1:n,1:m) + alpha_shift * ppsi(1:n,1:m)
END SUBROUTINE ch_psi_rest
```

### 6.2 Driver per channel

```fortran
SUBROUTINE solve_sternheimer_k(ik, nrhs, rhs, chi, conv)
  !! Solve Q(ik)(omega0 - H0(ik))Q(ik) chi = rhs  for nrhs right-hand sides.
  USE lr_module_shim, ONLY : ccgsolve_all          ! from LR_Modules/liblrmod.a
  IMPLICIT NONE
  INTEGER,     INTENT(IN)  :: ik, nrhs
  COMPLEX(dp), INTENT(IN)  :: rhs(:,:)
  COMPLEX(dp), INTENT(OUT) :: chi(:,:)
  LOGICAL,     INTENT(OUT) :: conv
  REAL(dp)    :: e0, h_diag(npwx*npol, nrhs), anorm
  COMPLEX(dp) :: freq_c
  INTEGER     :: kter
  ! ccgsolve_all(ch_psi, ccg_psi, e, d0psi, dpsi, h_diag, ndmx, ndim, ethr,
  !              ik, kter, conv_root, anorm, nbnd, npol, freq_c)  solves (H - e + Q + freq_c) x = b
  e0     = omega0                       ! REAL reference energy
  freq_c = CMPLX(0._dp, eta, dp)        ! complex shift: omega0 -> omega0 + i*eta  (research.md §12)
  CALL g2_kin(ik)                       ! set g2kin for h_psi / preconditioner at ik
  CALL build_h_diag(ik, e0, freq_c, h_diag)         ! precond ~ 1/(g^2 - e0 - freq_c)
  chi = (0._dp,0._dp)
  ! d0psi = -rhs  (since Q(omega0-H0)Q = -(H0-omega0) on Q)
  CALL ccgsolve_all(ch_psi_rest, ccg_psi, e0, -rhs, chi, h_diag, &
                    npwx*npol, ngk(ik)*npol, sternheimer_thr, ik, kter, conv, anorm, &
                    nrhs, npol, freq_c)
  ! enforce P|chi>=0 exactly on exit
  CALL apply_Qproj(ngk(ik), npol, ik, n_act(ik), evc_act(:,:,ik), chi)   ! per rhs
END SUBROUTINE solve_sternheimer_k
```

This is the **only** linear-algebra novelty; everything is per-$k$, well-conditioned
(`research.md` §7.1), and embarrassingly parallel over the $k$-grid (§11).

---

## 7. Stage D — $k$-decoupled rest dressing ladder (NEW, optional)

Bare ($V_{QQ}=0$) — i.e. **Stage C, the $m=0$ rung** — gives the coupling-second-order $\tilde V^{(2)}$.
To reach the *exact* dressed $\tilde V$ while staying $k$-decoupled, Neumann-iterate (`research.md` §9.1):
$$|\chi^{(0)}\rangle=G^R V_{QP}|n\rangle,\qquad |\chi^{(m)}\rangle=G^R\,V_{QQ}\,|\chi^{(m-1)}\rangle,\qquad
|\chi\rangle=\textstyle\sum_m|\chi^{(m)}\rangle .$$
Each step is **one $V_{QQ}$ matvec** (reuse §5's ΔV-action restricted to $Q$, *no $k$-coupling in the
solve*) followed by the **same** per-$k$ bare solve $A_0=Q(\omega_0-H_0)Q$.

```fortran
chi = chi0                                   ! = G^R Q ΔV |n>  (Stage C with bare A0)
DO mstep = 1, max_dress
   CALL apply_dV_to_ket(chi, vqq_chi)        ! V_QQ * chi  (ΔV-action, then apply_Qproj)
   DO ik = 1, nks
      CALL solve_sternheimer_k(ik, nrhs, vqq_chi(:,:,ik), dchi(:,:,ik), conv)  ! A0^{-1} on Q
   END DO
   chi = chi + dchi
   ratio = norm(dchi) / norm_prev            ! successive-ratio diagnostic (research.md §10)
   IF (ratio < dress_tol) EXIT               ! geometric rate ρ ~ ||V_QQ||/Δ
   norm_prev = norm(dchi)
END DO
```

(GMRES with $A_0$ as preconditioner and $V_{QQ}$ as the FFT matvec is the accelerated cousin; drop-in
once Neumann works.)

---

## 8. Stage E — assemble $\tilde V$ (Hermitian bilinear form)

Per $(n k_i, m k_f)$ active pair (`research.md` §7.3, bare: $V_{QQ}\to0$):
$$\tilde V_{m k_f,\,n k_i}=\underbrace{\langle m k_f|\Delta V|n k_i\rangle}_{=\,M\ (\text{EDI})}
\;+\;\big\langle\chi_{m k_f}\big|(\omega_0-H_0^Q-V_{QQ})\big|\chi_{n k_i}\big\rangle .$$

```fortran
! M-term: reuse EDI's coarse matrix element (ed_coarse_full_q kernel) -> edmat_bloch(m,n,ki,kf)
! rest-term: contract the stored responses over channels k'
DO ik_f = ...; DO ik_i = ...
  DO m = 1, n_act(ik_f); DO n = 1, n_act(ik_i)
     rest = (0._dp,0._dp)
     DO ikp = 1, nks                                   ! channels k'
        ! <chi_m| (omega0 - H0^Q - V_QQ) |chi_n>  at channel k'
        CALL ch_rest_metric(ikp, chi(:,m,ik_f,ikp), chi(:,n,ik_i,ikp), val)
        rest = rest + val
     END DO
     Vtilde_B(m,n,ik_i,ik_f) = edmat_bloch(m,n,ik_i,ik_f) + rest
  END DO; END DO
END DO; END DO
! Hermitize in (m k_f) <-> (n k_i): Vtilde <- (Vtilde + Vtilde^H)/2   (real omega0)
```

`ch_rest_metric` reuses `ch_psi_rest` (it *is* $\omega_0-H_0$ on $Q$, plus the $V_{QQ}$ term via §5's
matvec). At the bare level there is a cheaper identity (`research.md` §7.3 derivation):
$\langle\chi_m^{(0)}|(\omega_0-H_0^Q)|\chi_n^{(0)}\rangle=\langle\chi_m^{(0)}|s_n\rangle$, i.e. just overlap the
response with the source — no extra operator application.

---

## 9. Stage F — active-layer resummation $T_{PP}(\omega)$ (NEW, small dense solve)

Active space $A=\{(a,k):a\in\text{active}(k)\}$, dimension $N_A=\sum_k n_{\rm act}(k)$ (coarse grid:
$\sim$ bands$\times N_k$, e.g. $11\times144=1584$). $G^A(\omega)$ is **diagonal** in this basis from the
interpolated active eigenvalues; $\tilde V$ is the dense block from Stage E. *(This is the small **active**
manifold only — distinct from §5's dense **rest** BZ grid: the rest sum is the expensive BZ integral, the
active inversion is the cheap small block.)*

```fortran
SUBROUTINE active_tmatrix(omega, eta_a, Vtilde, eig_act, Tpp)
  !! Tpp = (I - Vtilde * G^A(omega))^{-1} * Vtilde   ,  G^A diagonal: 1/(omega - eig_act + i eta_a)
  USE kinds, ONLY : dp
  IMPLICIT NONE
  REAL(dp),    INTENT(IN)  :: omega, eta_a, eig_act(:)        ! eig_act(NA)
  COMPLEX(dp), INTENT(IN)  :: Vtilde(:,:)                     ! (NA,NA)  active basis
  COMPLEX(dp), INTENT(OUT) :: Tpp(:,:)
  COMPLEX(dp), ALLOCATABLE :: A(:,:), ga(:)
  INTEGER,     ALLOCATABLE :: ipiv(:)
  INTEGER :: NA, j, info
  NA = SIZE(eig_act)
  ALLOCATE(A(NA,NA), ga(NA), ipiv(NA))
  DO j = 1, NA
     ga(j) = 1._dp / CMPLX(omega - eig_act(j), eta_a, dp)
  END DO
  A = -MATMUL(Vtilde, DIAG(ga));  DO j=1,NA; A(j,j)=A(j,j)+1._dp; END DO   ! A = I - Vtilde G^A
  Tpp = Vtilde
  CALL ZGESV(NA, NA, A, NA, ipiv, Tpp, NA, info)              ! Tpp = A^{-1} Vtilde
END SUBROUTINE active_tmatrix
```

**Resolution choice (the main Layer-2 approximation — flag for validation):**
*baseline* = do this inversion on the **coarse** active space, then Wannierize $T_{PP}(k_i,k_f)$ and
interpolate to the fine grid (mirrors EDI's M-flow). *Advanced* = keep the static $\tilde V$ Wannier-
interpolated and solve the Dyson equation $T_{PP}=\tilde V+\tilde V G^A(\omega)T_{PP}$ iteratively on the
fine grid (matvec with interpolated $\tilde V$; never form the dense inverse), exploiting that $G^A(\omega)$
is dominated by the energy shell near $\omega$. Start with baseline; validate against advanced on a
medium grid.

---

## 10. Stage G — Wannierize, interpolate, transport (reuse EDI)

$\tilde V$ (and $T_{PP}$) carry the **same two-momentum, active-band index structure** as EDI's $M$, so the
Wannier interpolation is identical — but it is worth deriving the transform once to fix the gauge and
phase conventions and to see *why* it is the right object to interpolate.

### 10.1 Setup and the Wannier functions

Active Bloch states $|\psi_{n k}\rangle$ ($n\in A(k)$), eigenvalues $\varepsilon_{nk}$, on the coarse grid
$\{k\}$ ($N_k$ points). The disentanglement+localization rotation $U(k)$ (EDI's `cu`=`u_kc` from `filukk`)
defines the Wannier functions (EPW/W90 convention)
$$
|w_{iR}\rangle=\frac{1}{N_k}\sum_{k}\sum_{n\in A(k)} e^{-ik\cdot R}\,U_{ni}(k)\,|\psi_{nk}\rangle,
\qquad i=1\ldots N_W .
$$
The object produced by Layer 1 is the **effective potential in the active-Bloch basis**,
$\tilde V_{nm}(k_i,k_f)\equiv\langle\psi_{n k_i}|\tilde V|\psi_{m k_f}\rangle$ (Stage E), the downfolded
replacement for EDI's bare $M_{nm}(k_i,k_f)=\langle\psi_{n k_i}|\Delta V|\psi_{m k_f}\rangle$.

### 10.2 Bloch → Wannier (the double Fourier transform) — derivation

The effective potential between two Wannier functions is, inserting the definition above and its conjugate
$\langle w_{iR_e}|=\frac1{N_k}\sum_{k_i,n}e^{+ik_i\cdot R_e}U^{*}_{ni}(k_i)\langle\psi_{n k_i}|$,
$$
\tilde V_{ij}(R_e,R_p)\equiv\langle w_{iR_e}|\,\tilde V\,|w_{jR_p}\rangle
=\frac{1}{N_k^{2}}\sum_{k_i,k_f}\sum_{n,m}
e^{+ik_i\cdot R_e}\,e^{-ik_f\cdot R_p}\;
U^{*}_{ni}(k_i)\,\tilde V_{nm}(k_i,k_f)\,U_{mj}(k_f),
$$
i.e., in matrix form (Wannier indices $i,j$; band indices contracted),
$$
\boxed{\;\tilde V(R_e,R_p)=\frac{1}{N_k^{2}}\sum_{k_i,k_f}
e^{+ik_i\cdot R_e}\,e^{-ik_f\cdot R_p}\;
U^{\dagger}(k_i)\,\tilde V(k_i,k_f)\,U(k_f)\;}
\tag{G.1}
$$
This is **exactly EDI's `edbloch2wane` double FT** (phases $+k_i\!\cdot\!R_e$, $-k_f\!\cdot\!R_p$; normalization
$1/N_k^2$; Wannier-gauge rotation $U^\dagger(k_i)(\cdot)U(k_f)$) — only the matrix fed in changes
$M\!\to\!\tilde V$. The very same lines apply to $T_{PP}$: $T_{PP}(R_e,R_p)$ is (G.1) with
$\tilde V(k_i,k_f)\!\to\!T_{PP}(k_i,k_f)$.

### 10.3 Wannier → Bloch (interpolation to the fine grid) — derivation

Inverting (G.1) on the Wigner–Seitz supercell with degeneracy weights $N_{R}$ gives, for *arbitrary*
fine-grid momenta,
$$
\boxed{\;\tilde V^{W}(k_i,k_f)=\sum_{R_e,R_p}
\frac{e^{-ik_i\cdot R_e}}{N_{R_e}}\,
\frac{e^{+ik_f\cdot R_p}}{N_{R_p}}\;
\tilde V(R_e,R_p)\;}
\tag{G.2}
$$
(EDI's `edmatwan2bloch_2d`: conjugate phase on the bra $k_i$, direct phase on the ket $k_f$, divided by the
WS degeneracies). (G.2) returns the object in the **Wannier gauge**; rotate to the band gauge at the
interpolated $k$ with the Hamiltonian eigenvectors $U^{H}(k)$ from diagonalizing the interpolated
$H_W(k)$ (`hamwan2bloch_with_evec`):
$$
\tilde V_B(k_i,k_f)=U^{H\dagger}(k_i)\,\tilde V^{W}(k_i,k_f)\,U^{H}(k_f).
\tag{G.3}
$$
**Locality / why interpolation is controlled.** $\tilde V=V_{PP}+V_{PQ}\mathcal G^Q V_{QP}$ is the bare local
defect block plus a *localized* rest self-energy (the rest response $|\delta\psi^R\rangle$ is short-ranged
for a gapped rest), so $\tilde V(R_e,R_p)$ decays in $|R_e|,|R_p|,|R_e-R_p|$ just like EDI's $M(R,R')$
(`diagonal_approx.md`). Monitor a `decay.V` plot (analog of EDI's `decay.M` from `edbloch2wane`).

### 10.4 Two routes to $T_{PP}$ on the fine grid

| | route | when |
|---|---|---|
| **B1** | **Wannierize $T_{PP}$ directly**: build $T_{PP}(k_i,k_f)$ on the *coarse* grid (active inversion §9), then apply (G.1)–(G.3) with $\tilde V\!\to\!T_{PP}$. Identical code path. | default *iff* $T_{PP}(R_e,R_p)$ decays |
| **B2** | **Wannierize $\tilde V$ (short-ranged), resum on the fine grid**: interpolate $\tilde V$ via (G.2)–(G.3), then solve the active Dyson equation (below) iteratively per energy. | when $T_{PP}$ is long-ranged (near-resonant / bound states) |

**Caveat for B1.** The resummation $T_{PP}=[1-\tilde V G^A]^{-1}\tilde V$ can be **longer-ranged** than
$\tilde V$: multiple scattering and any near-shell resonance spread the effective range in $(R_e,R_p)$.
So **check the `decay.T` of $T_{PP}(R_e,R_p)$** before trusting B1; if it does not decay within the coarse
WS cell, use B2 (which only ever interpolates the short-ranged $\tilde V$).

**B2 Dyson equation (active space, fine grid).** With $G^A(\omega)$ diagonal in the band/$k$ basis,
$$
T_{PP}(k_i,k_f;\omega)=\tilde V_B(k_i,k_f)
+\frac{1}{N_k^{\rm f}}\sum_{k'}\sum_{a\in A(k')}
\tilde V_B(k_i,k')\,\frac{|a k'\rangle\langle a k'|}{\omega-\varepsilon_{a k'}+i\eta}\,T_{PP}(k',k_f;\omega),
\tag{G.4}
$$
solved by fixed-point/GMRES iteration using interpolated $\tilde V_B$ as a matvec — never forming the dense
inverse. The $1/N_k^{\rm f}$ turns the $k'$-sum into a BZ average; the resolvent is dominated by the
energy shell $\varepsilon_{ak'}\approx\omega$.

### 10.5 Transport (reuse EDI unchanged)

The on-shell $T$-matrix replaces $M$ in Fermi's golden rule,
$$
\frac{1}{\tau_{nk}}=\frac{2\pi}{\hbar}\,n_d\,\frac{1}{N_k}\sum_{m,k'}
\big|\,T_{PP,nm}(k,k';\,\omega{=}\varepsilon_{nk})\,\big|^{2}\,
\delta(\varepsilon_{nk}-\varepsilon_{mk'}),
$$
i.e. in `compute_transport` the only change is `ABS(edmatf_b)**2 → ABS(Tpp_b)**2`; IBZ symmetry,
SERTA/MRTA, `delta_weights`, and the Fermi-level bisection are reused verbatim. The Born limit
$T_{PP}\to\tilde V\to M$ reproduces today's EDI mobility (test **T1**).

```fortran
! ---- coarse: Layer 1 (+ optional Layer 2 for route B1) ----
CALL edbloch2wane(nbnd_kept, nbndsub, nks, nkstot, xk, cu, cuq, &   ! (G.1) double FT  [reuse EDI]
                  Vtilde_B, nrr, irvec, wslen, Vtilde_W)            ! -> Vtilde(R_e,R_p); writes decay.V
! route B1: replace Vtilde_B by Tpp_B (active_tmatrix on coarse grid) before the FT, check decay.T

! ---- fine grid: per (k_i,k_f) ----
CALL get_cfac(nrr, irvec, xk_i_cryst, cfac_ki)                      ! exp(i k_i·R)   [reuse EDI]
CALL get_cfac(nrr, irvec, xk_f_cryst, cfac_kf)
CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, Vtilde_W, cfac_ki, cfac_kf, Vtilde_Wf)  ! (G.2) [reuse]
! band gauge (G.3): Vtilde_Bf = U^H(k_i)^dagger · Vtilde_Wf · U^H(k_f)   (evec from hamwan2bloch_with_evec)
!   route B1: same for Tpp_W -> Tpp_Bf
!   route B2: solve (G.4) for Tpp_Bf using Vtilde_Bf as the matvec

! ---- golden rule (reuse transport.f90 kernel) ----
!   inv_tau += twopi * n_d * wqf * ABS(Tpp_Bf(n,m))**2 * w_delta      (was ABS(edmatf_b)**2)
```

---

## 11. Parallelization (mirror EDI + the codex run)

- **Outer:** distribute *source states* $(\text{initial }k_i,\ \text{active band }n)$ across MPI images /
  pools (EDI's `mp_pools` + panel-broadcast pattern). The reference codex run did fixed $k_i$, rest/final
  $k=1\!:\!144$, **1584 RHS solves** in $\sim$2h27m on 72 image-MPI ranks — EDT reproduces this layout.
- **Inner:** per source, the rest-$k'$ channel solves (§6, over the full-BZ rest grid `rest_nk*`) are
  independent → vectorize as `nrhs` columns to `ccgsolve_all`, or split across the pool. The rest grid is
  typically the dominant cost (sources × rest-$k'$ × solve), so its density is the main convergence/cost knob.
- **Reductions:** the rest-term contraction (§8) sums over channels with `mp_sum(inter_pool_comm)`, like
  `edbloch2wane`.
- Cache per-$k$ real-space wavefunctions and `becd`/`becp` once (EDI already does this in
  `ed_coarse_full_q`).

---

## 12. Input namelist `&edt_nml` (extends EDI's `&edinput_nml`)

| Key | Type | Default | Meaning |
|---|---|---|---|
| *(all EDI keys)* | — | — | `edi_prefix/outdir`, `potfile_*`, `pot_align`, `coarse_nk*`, `fine_nk*`, Wannier keys, transport keys, `range_sep`, … |
| `do_tmatrix` | logical | `.false.` | master switch for the EDT path |
| `active_win_min/max` | real (eV) | — | active window $[\,]$ defining $A(k)$ (defaults to transport window) |
| `rest_nk1/2/3` | integer | = coarse | **rest BZ grid** for $G^R$ / the Sternheimer sum — spans the *full* BZ, $\ge$ coarse, densified to convergence (§5 callout). Not the supercell folding set. |
| `omega0` | real (eV) | window center | Sternheimer reference energy |
| `eta` | real (eV) | `0.01` | imaginary shift for `ccgsolve_all` (and metallic case) |
| `sternheimer_thr` | real | `1e-10` | linear-solver residual (codex reached `9.98e-11`) |
| `dress_order` | integer | `0` | $V_{QQ}$ ladder steps; `0` = coupling-2nd-order $\tilde V^{(2)}$ |
| `dress_tol` | real | `1e-3` | successive-ratio stop for the ladder (§7) |
| `active_resum` | logical | `.true.` | do Layer 2 ($T_{PP}$); `.false.` ⇒ stop at $\tilde V$ |
| `resum_grid` | char | `'coarse'` | `'coarse'` (invert+interpolate) or `'fine'` (iterative Dyson) |
| `omega_grid_*` | real | transport win | frequencies for $T_{PP}(\omega)$ |
| `rest_split` | char | `'complex'` | `'complex'` (ccgsolve, $\omega_0+i\eta$) or `'pm'` (split $R^\pm$, cgsolve) |

---

## 13. Build system

`edt/src/makefile` mirrors EDI's but additionally links `LR_Modules` (for the Sternheimer solver) and the
EDI objects:

```make
include ../../make.inc
QEMODS = ../../Modules/libqemod.a ../../upflib/libupf.a \
         ../../KS_Solvers/libks_solvers.a ../../LR_Modules/liblrmod.a \   # <-- NEW: cgsolve/ccgsolve
         ../../FFTXlib/src/libqefft.a ../../UtilXlib/libutil.a \
         ../../XClib/xc_lib.a ../../LAXlib/libqela.a
PWOBJS = ../../PW/src/libpw.a
EDIOBJS = ../../edi-code/src/ed_coarse.o ../../edi-code/src/range_sep.o \
          ../../edi-code/src/get_betavkb.o ../../edi-code/src/wan2bloch.o \
          ../../edi-code/src/edbloch2wan.o ../../edi-code/src/transport.o \
          ../../edi-code/src/edic_mod.o ../../edi-code/src/wann_common.o   # reuse compiled EDI
edt.x : $(EDT_OBJS) ; $(LD) -o $@ $(EDT_OBJS) $(EDIOBJS) $(PWOBJS) $(QEMODS) $(QELIBS) $(LIBS)
```

(If EDI's `.o` ABI drifts, vendor the handful of reused EDI source files into `edt/src/` instead.)

---

## 14. Validation & tests (feed the website Test Catalog)

| # | Test | Pass criterion | Maps to |
|---|---|---|---|
| T1 | **Born limit** | `dress_order=0, active_resum=.false.` and rest-term zeroed ⇒ $\tilde V \equiv M$ to $\lesssim10^{-12}$ Ry vs `ed_coarse_full_q` | regression vs EDI |
| T2 | **Sternheimer vs explicit sum** | $\tilde V^{(2)}$ from `ccgsolve` matches the explicit finite-rest-band sum $\sum_{r\in R}V_{Pr}V_{rP}/(\omega_0-\varepsilon_r)$ as band cutoff → all | `research.md` §3.2; codex "explicit vs QE Sternheimer" |
| T3 | **Solver health** | all RHS converge below `sternheimer_thr`; report residual, $A_0(k)$ smallest-$|{\rm eig}|$ / condition number vs $\lVert V\rVert$ | `research.md` §10 |
| T4 | **Ladder convergence** | successive ratio $\to\rho\sim\lVert V_{QQ}\rVert/\Delta\ll1$; $\tilde V$(dressed) stable vs `dress_order` | §7, §10 |
| T5 | **Gauge sanity** | $\lVert U^\dagger U-I\rVert\!\lesssim\!10^{-13}$; $\tilde V$ invariant under $U^\dagger(\!\cdot\!)U$ | EDI `filukk` |
| T6 | **$k$-MPI invariance** | result independent of pool/image count to machine precision | codex image-MPI check |
| T7 | **Active resummation** | coarse-grid $T_{PP}$ (B1) vs fine-grid iterative Dyson (B2) agree within tol; $T_{PP}\to\tilde V$ as $G^A\to0$; `decay.V`/`decay.T` decide B1-vs-B2 | §9, §10.4 |
| T9 | **Rest BZ-grid convergence** | $\lVert\tilde V\rVert$ (and mobility) converge as `rest_nk*` densifies over the *full* BZ; result is **independent of the supercell folding set** | `research.md` §2; §5 callout; codex rest-k table |
| T8 | **Transport** | mobility with $|T_{PP}|^2$ vs EDI $|M|^2$; quantify beyond-Born shift | `transport.f90` |

Each test emits an `ionode` report + a small CSV (numbers only), summarized as a row in the
`claude-sternheimer` site Test Catalog (see repo `CLAUDE.md` for the page recipe). **Never** publish
wavefunctions/cubes/`.npy`/logs.

---

## 15. Milestones

- **P0 — scaffold:** `edt.f90` + `&edt_nml`; call EDI Stage A; define $A(k)$, `apply_Qproj`; print active/rest audit (counts, $\omega_0$, gap $\Delta$). *(no solve yet)*
- **P1 — source:** `edt_source` (`dV_local_ket`+`dV_nonlocal_ket`); **T1** Born-limit check that $\langle m k_f|$source$\rangle$ reproduces $M$.
- **P2 — Sternheimer:** `ch_psi_rest` + `solve_sternheimer_k` (ccgsolve); single $(k_i,n)$ smoke; **T2/T3**.
- **P3 — $\tilde V$:** Stage E bilinear assembly + Hermitize; **T5**; full $\tilde V$ for fixed $k_i$ over the rest grid (the codex "1584 RHS" milestone), then **T9** rest-BZ-grid convergence.
- **P4 — dressing:** `edt_dress` ladder; **T4**.
- **P5 — active layer:** `active_tmatrix`; **T7**.
- **P6 — interpolate+transport:** reuse `edbloch2wane`/`edmatwan2bloch_2d`/`compute_transport`; **T6/T8**; first beyond-Born mobility.

## 16. File checklist

**New:** `edt.f90`, `edt_input.f90`, `edt_partition.f90`, `edt_source.f90`,
`edt_sternheimer.f90`, `edt_dress.f90`, `edt_vtilde.f90`, `edt_active.f90`, `src/makefile`,
`tests/` (T1–T8 inputs + report scripts).
**Reused from EDI (link or vendor):** `ed_coarse.f90`, `range_sep.f90`, `get_betavkb.f90`,
`wan2bloch.f90`, `edbloch2wan.f90`, `transport.f90`, `delta_weights.f90`, `bz_symmetry.f90`,
`edic_mod.f90`, `wann_common.f90`, `global_var.f90`, `edi_input.f90` (extended).
**QE libs:** `libpw.a`, `libqemod.a`, `liblrmod.a` (**new link**), `libks_solvers.a`, `libqefft.a`, `libupf.a`, `xc_lib.a`, `libutil.a`, `libqela.a`.

---

## 17. Detailed TODO checklist (all phases & tasks)

Granular, checkable expansion of the §15 milestones. **Legend:** `[ ]` to do · `[x]` done · **(§n)**
stage/section · **(Tn)** validation test (§14) · `[file]` target source. Phases are ordered for
incremental validation — **each phase ends in a test gate that must pass before the next**. Nothing here
is implemented yet (this is the plan only).

### P0 — Scaffold, build system, inputs & active/rest partition  ✅ **DONE (2026-05-30)**
*Validated on MoS₂ (150-band NSCF, 12×12): active = 11 Wannier valence bands, rest = 133 per k
(19152 total); Wannier interp vs NSCF = 1.8×10⁻⁵ eV; ΔV cube load (S-vacancy, 240×240×300) OK;
MPI pool-invariant (npool=2 ≡ npool=1, early T6 check).*
- [x] Decide packaging: **sibling plug-in `edt/src/` in the repo, linking compiled EDI objects + QE libs by absolute path** (§2, §13)
- [x] `edt/src/` tree; `src/makefile` links `libpw.a`, `libqemod.a`, **`liblrmod.a`** (new), `libks_solvers.a`, FFT/upf/xc/util/lax via `--start-group` (§13) `[makefile]`
- [x] Link resolves `h_psi`, `ccgsolve_all`, `cgsolve_all`, `cg_psi` from `LR_Modules` (edt.x links clean; solver *called* in P2) (§6)
- [x] `&edt_nml` with all keys + `mp_bcast` every key (§12) `[edt_input.f90]`
- [x] Main program: read+broadcast input, run Stage A + audit `[edt.f90]`
- [x] Stage A loaders (reuse EDI): `read_file`, `load_pot_from_file`, `build_vcolin_aligned`/`_corealign`, `read_hr_file`, Wannier-interp sanity (1.8e-5 eV) (§4). **NOTE:** `read_filukk_edi` overruns when NSCF nbnd ≠ Wannierization nbnd → wrote minimal robust reader `edt_wannier.f90::edt_read_filukk` (`nbndep`/`ibndkept`/`u_kc` only)
- [x] Optional range separation: `compute_range_separation` wired + guarded (coded; not exercised — no `rhofile_d/p` in the MoS₂ test) (§4; `research.md` §16.5)
- [x] Build active set `is_active(ibnd,ik)`, `n_act(ik)` from window; set $\omega_0$ (=VBM default), compute gap $\Delta$ `[edt_partition.f90]`
- [x] `apply_Qproj` — Gram–Schmidt against active states (NC: $S=1$) (§3) `[edt_partition.f90]` (defined; first *exercised* in P1)
- [x] Active/rest audit print (ionode): $N_A$, $N_R$, $\omega_0$, $\Delta$, window, per-$k$ counts
- [x] **Gate:** `edt.x` builds and runs the audit on MoS₂ (no solve yet) ✓

### P1 — Stage B: Sternheimer source $|s\rangle=Q\,\Delta V\,|nk_i\rangle$  ✅ **DONE (2026-05-30)**
*Full source ket (local + nonlocal KB, defect−pristine) reproduces EDI's matrix-element kernels to
**max|ΔM| = 2.3×10⁻¹³ Ry** (MoS₂, q=(0,1/12,0), 5 bands; |M_loc|≈0.65, |M_nl|≈0.12 Ry; nkb_d=1926,
nkb_p=1944 — the S-vacancy removes exactly 18 KB projectors).*
- [x] `build_V_folded` (exact $q$, no BZ wrap) + local ΔV ket via `invfft`/`fwfft` on the $k'$ G-set (§5.1) `[edt_source.f90]` — T1-local 2.2e-13
- [x] nonlocal ket: `get_betavkb`($k_a$,$k_b$) + $\langle\beta|\psi\rangle$ + `dvan` contraction, defect − pristine (§5.2) `[edt_source.f90]` — T1-nonlocal 1.0e-13
- [x] **(T1) Gate — Born limit:** full local+nonlocal source-ket vs EDI kernels = **2.3e-13** ✓ (`test_source`; reuses EDI's `get_betavkb`/`dvan` verbatim, so it validates against EDI's exact kernels)
- [ ] *(folds into P2)* refactor `test_source` → reusable `dV_*_ket` + per-channel **full-BZ assembly** of $|s\rangle$ with `apply_Qproj` — built where the RHS is consumed (P2)
- [ ] *(deferred opt.)* `gk_sort` for fine rest-$k'$ (`rest_nk*`≠coarse) + ψ/`becp` caching (§11)
- [ ] *(nice-to-have)* end-to-end cross-check of full $M$ vs an actual `ed_coarse_full_q` run (current check is kernel-level)

### P2 — Stage C: bare per-$k'$ Sternheimer solve  ✅ **DONE (2026-05-30)**
*Self-written projected (Jacobi-PCG) solve of $[(H_0-\omega_0)+\alpha P_{\rm act}]\chi=s$ using QE `h_psi`
as matvec (rest is above ω₀ here ⇒ positive-definite). MoS₂ (isrc=17, ω₀=−1.09 eV):*
| | CG iters / resid | explicit cutoff 30→150 | **Sternheimer (all bands)** |
|---|---|---|---|
| q=0 | 60 / 3.8e-11 | −0.095→−0.233 Ry | **−0.301 Ry** |
| q≠0 | 76 / 6.5e-11 | −0.505→−0.668 Ry | **−0.733 Ry** |
*The explicit band sum monotonically **approaches** the Sternheimer value; Sternheimer captures the full
all-band result that explicit-150 misses by ~9–22% (the high-band tail) — the payoff of the method.*
- [x] **h_psi gate:** $\langle\psi_{nk}|H_0|\psi_{nk}\rangle=\varepsilon_{nk}$ to **8.9×10⁻¹⁰ eV** (`edt_set_vrs`+`hpsi_setup_k`+`test_hpsi_eigen`; validates vrs/vkb/g2kin) `[edt_sternheimer.f90]`
- [x] `solve_rest_cg` — projected Jacobi-PCG, `h_psi` matvec + $\alpha P_{\rm act}$ shift; `build_source_ket`+`apply_Qproj` make the RHS (§6) `[edt_sternheimer.f90]`
- [x] `rest_channel_compare` — explicit (`explicit_rest_channel`) vs Sternheimer at a channel `[edt_sternheimer.f90]`
- [x] **(T2)** explicit sum → Sternheimer all-band value ✓ (monotone approach + correct tail)
- [x] **(T3)** solver health: CG converges, residual ~$10^{-11}<$ `sternheimer_thr` ✓
- [ ] *(deferred):* full QE `ccgsolve_all` complex-shift path for the **indefinite/metallic** case (§12) — here the operator is definite so PCG suffices; the `ch_psi`/`euc` interface is mapped out.
- [x] **k′-sum normalization resolved:** the $k'$-sum is the **BZ-integral discretization** ⇒ carries **$1/N_k$** ($\mathbb 1=\frac1{N_k}\sum_{nk}|nk\rangle\langle nk|$, QE per-cell states), so $\Sigma_{mn}=\frac1{N_k}\sum_{k'}\sum_r MM/(\omega_0-\varepsilon)$. This turns the naive $\sum_{k'}\!\approx\!-70$ Ry into the physical $\Sigma_{nn}\!\approx\!-0.5$ Ry (≈ Born $M$). $G^A$ also carries $1/N_k$; ties to the §5 rest-grid BZ convergence and EDI's golden-rule $1/N_k$. See the [Implementation Note](note-kprime-normalization.html).
- [x] **no $N_{\rm sc}$ factor:** for converged $\Delta V$, $M$ (hence $\tilde V$, $\Sigma$) is **supercell-independent** — EDI's primitive-cell measure $\frac1{N_{\rm nnr}}$ on the *localized* folded $\Delta V$ only sees the defect region (concentration handled by $n_d$, not $M$). Two independent convergence axes: supercell size ($\Delta V$) and rest BZ grid $N_k$ ($\frac1{N_k}\sum\!\to\!\int$).
- [ ] **(→ P3) confirm:** per-channel closure $\sum_{n'}|M|^2\!=\!\lVert s\rVert^2$; full $\frac1{N_k}\sum_{k'}\sum_{n'}|M|^2\!=\!\langle\Delta V^2\rangle$ (both supercell-independent); $M$ invariant across supercell sizes; Born-limit mobility == EDI (gold-standard anchor).

### P3 — Stage E: assemble the coupling-second-order $\tilde V^{(2)}$  ✅ **done (diagonal + full block)**
*`vtilde_diag_full` assembles the **diagonal** $\tilde V_{nn}=M_{nn}+\Sigma_{nn}$ over the full BZ with the
$1/N_k$ measure, and runs the closure check. MoS₂ (isrc=17 VBM, ki=1, $N_k$=144):*
| $M_{nn}$ (Born) | $\Sigma_{nn}$ explicit-150 | $\Sigma_{nn}$ **Sternheimer** | **$\tilde V_{nn}$** | $\langle\Delta V^2\rangle$ (norm / 150-band) |
|---|---|---|---|---|
| +0.701 Ry | −0.773 Ry | **−0.819 Ry** | **−0.117 Ry** | 1.756 / 1.529 Ry² (ratio 0.871) |
*$1/N_k$ **confirmed**: $\langle\Delta V^2\rangle\!\sim\!\mathcal O(1)$ Ry² (not $144\times$); band completeness 0.871
(deficit = high-band tail, recovered by Sternheimer). $|\Sigma_{nn}|>M_{nn}$, sign-flipping ⇒ strong
beyond-Born regime (Born fails, T-matrix essential).*
- [x] rest-term contracted over the full-BZ channels with $1/N_k$ (BZ measure); $\Sigma_{nn}=\frac1{N_k}\sum_{k'}(-\langle s|\chi\rangle)$ `[edt_sternheimer.f90::vtilde_diag_full]`
- [x] Born $M_{nn}$ term (q=0 channel of the source ket); $\tilde V_{nn}=M_{nn}+\Sigma_{nn}$
- [x] **closure check** confirms $1/N_k$: $\frac1{N_k}\sum_{k'}\sum_{n'}|M|^2=\langle\Delta V^2\rangle\approx1.76$ Ry² ($\mathcal O(1)$), band completeness 0.871
- [x] **full off-diagonal block** $\tilde V_{(mk_f),(nk_i)}$ (all 1584 active states) + Hermitize `[edt_sternheimer.f90::vtilde_block_mpi]` — **pool-parallel** over the $k'$-sum (`-nk N`: 1 rank/pool ⇒ full G-grid, native `h_psi` at pool-local $k'$; gather $\varepsilon$/active wfc to all ranks; Born ZGEMM + Q-proj + Sternheimer + $\Sigma$ ZGEMV; `mp_sum` over `inter_pool_comm`; write `vtilde_block.dat`). MoS₂ **1584×1584** on 1 node / 36 ranks, 2 h 11 m.

  | per-rank $H_0$ gate | Hermiticity $\lVert\tilde V-\tilde V^\dagger\rVert$ | in-situ band-17/k=1 ($M$,$\Sigma$,$\tilde V$) |
  |---|---|---|
  | 6.0e-10 eV (all ranks) | **8.96e-12** (pre-symmetrization) | +0.70148, −0.81860, **−0.11712 Ry** — matches the diagonal run to 6 digits |

  *Confirms the cross-pool assembly is correct. Strong beyond-Born across the manifold (k=1: band 10 $\tilde V$=−0.234 Ry; band 7 −0.008 Ry from $M$=+0.49/$\Sigma$=−0.50 near-cancellation); degeneracies and BZ symmetry preserved. Memory-bound at ~3.4 GB/rank (supercell KB projectors $n_{kb}$=1926/1944), so ~40 ranks/257 GB node.*
- [x] **(T6)** MPI/pool correctness: full block is Hermitian to 9e-12 and its band-17/k=1 diagonal reproduces the single-rank value across 36 pools
- [ ] **(T5)** gauge: $\|U^\dagger U-I\|$; $\tilde V$ invariant under $U^\dagger(\cdot)U$
- [ ] **(T9)** rest BZ-grid convergence: $\|\tilde V\|$ vs `rest_nk*` densifying (full BZ)

### P4 — Stage D: $k$-decoupled $V_{QQ}$ dressing ladder (optional)
- [ ] Neumann ladder $\chi^{(m)}=G^R V_{QQ}\chi^{(m-1)}$; $V_{QQ}$ via the ΔV-action $\oplus\,Q$ (§7) `[edt_dress.f90]`
- [ ] Reuse `solve_sternheimer_k` (same $A_0$) per rung — only the RHS changes
- [ ] Successive-ratio diagnostic + `dress_tol` stop; honor `dress_order` ($0$ = keep $\tilde V^{(2)}$)
- [ ] (optional) GMRES variant: $A_0$ preconditioner, $V_{QQ}$ matvec
- [ ] **(T4) Gate — ladder convergence:** ratio $\to\rho\sim\|V_{QQ}\|/\Delta$; $\tilde V$ stable vs `dress_order`
- [ ] Deep/bound-level guard: detect $\rho\gtrsim1$ / $A_0$ near-singular → advise enlarging the active window (`research.md` §10)

### P5 — Stage F: active-space resummation $T_{PP}(\omega)$
- [ ] `active_tmatrix`: build diagonal $G^A(\omega)$, solve $[1-\tilde V G^A]^{-1}\tilde V$ via `ZGESV` (§9) `[edt_active.f90]`
- [ ] Pack the coarse active space (band$\times k$ indexing) for $\tilde V$ and `eig_act`
- [ ] $\omega$-grid loop (`omega_grid_*`); shell smearing `eta_a`
- [ ] **(T7a) Gate — sanity:** $T_{PP}\to\tilde V$ as $G^A\to0$

### P6 — Stage G: Wannierize, interpolate, transport
- [ ] Bloch→Wannier (G.1): reuse `edbloch2wane` on $\tilde V_B\to\tilde V(R_e,R_p)$; write `decay.V` (§10.2)
- [ ] Wannier→Bloch (G.2)+(G.3): reuse `edmatwan2bloch_2d` + `hamwan2bloch_with_evec` band gauge (§10.3)
- [ ] Route **B1**: Wannierize $T_{PP}$ directly — **check `decay.T`** before trusting (§10.4)
- [ ] Route **B2**: fine-grid Dyson (G.4) iterative solve (matvec with interpolated $\tilde V_B$) (§10.4)
- [ ] `resum_grid` switch (coarse B1 / fine B2)
- [ ] Transport hook: feed $|T_{PP}|^2$ into `compute_transport` (`ABS(edmatf_b)**2 → ABS(Tpp_b)**2`) (§10.5)
- [ ] Reuse `delta_weights`, `bz_symmetry`, Fermi-level bisection unchanged
- [ ] **(T7b)** B1 vs B2 agree within tol on a medium grid
- [ ] **(T8) Gate — transport:** mobility with $|T_{PP}|^2$ vs EDI $|M|^2$; quantify the beyond-Born shift

### P7 — Validation suite, long-range, performance, docs/release
- [ ] Assemble **T1–T9** into ionode reports + numbers-only CSVs (EDI diagnostic style)
- [ ] Add each test as a Test Catalog row on the site (repo `CLAUDE.md` recipe); flip badges `plan → ok/prod`
- [ ] Long-range: verify $\Delta V_{\rm SR}$ (Sternheimer) + analytic $M^{\rm LR}$ (`lr_matelem`) recombine (`research.md` §16.5)
- [ ] Performance: vectorize `nrhs` → `ccgsolve_all`; cache $\psi$/`becp`; profile rest-grid cost; MPI scaling (§11)
- [ ] Convergence study: `rest_nk*`, `dress_order`, `eta`, active-window sensitivity → short writeup
- [ ] End-to-end regression: the Born-limit pipeline reproduces EDI mobility
- [ ] (stretch) SOC path (`dvan_so`, `npol=2`) parity with EDI nonlocal
- [ ] (stretch) metallic host: complex $\omega_0+i\eta$ absorptive $\tilde V$ (`research.md` §12)

---

*Cross-references are to `research.md` (theory). This plan reuses the EDI implementation for the
difference potential, KB projectors, Wannier rotation/interpolation, and transport; the new code is the
rest-space Sternheimer solve, the $k$-decoupled $V_{QQ}$ ladder, the $\tilde V$ assembly, and the small
active-space resummation.*
