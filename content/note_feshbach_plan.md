# All-order Feshbach Sternheimer: implementation plan for the rest-space solve (v2)

The [ladder page](sternheimer-ladder.html) §5 measured the order-by-order expansion to **diverge** on the full rest space ($r_3=1.54$), and §6 derived the all-order Feshbach route that inverts directly instead of expanding. This page is its **code implementation plan (v2)**: taking a **static reference** for the rest contribution ($\omega\to\omega_0$, chosen inside the gap), the all-order $\Sigma$ becomes a single Sternheimer linear solve coupling all $\mathbf k$,

$$A\,|X_b\rangle=|S_b\rangle,\qquad A=Q\,(H_0+\Delta V-\omega_0)\,Q+\alpha P,\qquad
|S_b\rangle=Q\,\Delta V\,|b\rangle,\qquad \Sigma_{ab}(\omega_0)=-\langle S_a|X_b\rangle .$$

The only new piece of physics machinery is the **cross-channel $\Delta V$ operator** inside the matvec (`apply_dV`, a supercell implementation). The first version (v1) was written and run through **isolated acceptance — which caught a convention bug of order $1/N_k^2$**. This page starts with the post-mortem, then gives the corrected design: an acceptance ladder built on **element-by-element regression against the already-verified $M$ matrix**, a zero-communication source-parallel layout, an honest analysis of the operator's indefiniteness, and costs anchored to measurements.

## 1. Post-mortem of the v1 acceptance: three structural root causes

v1 (a pre-stored $\Delta V$ field plus a supercell FFT operator) compiled and ran cleanly, and isolated acceptance used `do_sigma3`'s folded $\Sigma^{(3)}$ as the reference:

| quantity | value |
|---|---|
| $\Sigma^{(3)}$, folded reference (raw) | $1.526\times10^{4}$ |
| $\Sigma^{(3)}$, v1 supercell FFT (raw) | $6.70\times10^{-1}$ |
| **ratio** | $4.4\times10^{-5}\approx0.91/N_k^2$ (against $1/N_k^2=4.82\times10^{-5}$, $N_k{=}144$) |

A clean $\sim1/N_k^2$ plus a ~9% residual ⇒ not a single missing factor, but a **structural error at the level of conventions**. Three root causes:

1. **The intra-cell Bloch phase $e^{i\mathbf k\cdot\mathbf r}$ was dropped.** v1's expansion/fold-back used $\sum_k u_k(\mathbf r)\,e^{i\mathbf k\cdot\mathbf R}$, but the real wavefunction is $\psi(\mathbf r{+}\mathbf R)=\sum_k e^{i\mathbf k\cdot(\mathbf r+\mathbf R)}u_k(\mathbf r)$ — **both** the intra-cell and inter-cell phases are needed. `build_V_folded`'s phase convention (`arg=irx*d1`, with positions including the intra-cell fraction) is exactly the complete $\mathbf q\cdot(\mathbf r{+}\mathbf R)$; v1 supplied only the inter-cell half, destroying coherence.
2. **The way the $\Delta V$ supercell field was assembled is structurally wrong.** v1 built it by $\mathbf q$-transforming the 144 folded potentials $V_f(\mathbf q)$ — but $V_f$'s phase contains a continuous $\mathbf r$, so the $\mathbf q$-sum is not a clean $\delta$ and produces aliasing (that ~9%). **The correct route never goes through $V_f(\mathbf q)$ at all**: the supercell field of $\Delta V$ is simply `V_colin` (the aligned cube data), already in memory, and needs only one grid rearrangement.
3. **The support region was wrong.** The support of $\Delta V$ is the **6×6 cube region (36 unit cells)**, not the 144 cells of the 12×12 BvK (the BvK cell = 2×2 defect supercells, and the potential is nonzero in only one of them). v1 spread the field over 144 cells — both wrong and 4× more expensive.

**Methodological lesson:** scalar-level acceptance (one $\Sigma^{(3)}$ number) can detect an error but cannot localize it. v2's foundation is instead **element-by-element** regression (V0 in §4).

## 2. The normalization contract (fixed once, and appearing only at file boundaries)

In the BvK-normalized basis $|kG\rangle$, a state is a vector of evc coefficients stacked over channels, and $\langle k'G'|\Delta V|kG\rangle$ **carries its own $1/N_k$**. The physical operator:

$$\big[\Delta V_{\rm phys}X\big]_{k'}=\frac1{N_k}\sum_{\mathbf k}V_f^{\,\mathbf k-\mathbf k'}\!\circ u_{\mathbf k},
\qquad V_f^{\mathbf q}(\mathbf r)=\sum_{\mathbf R\in\text{cube}}\Delta V(\mathbf r{+}\mathbf R)\,e^{i\mathbf q\cdot(\mathbf r+\mathbf R)} .$$

The consistency chain (all verifiable): $M_{\rm phys}=M_{\rm raw}/N_k$ (the project's existing contract, with $H_{\rm eff}$ already verified against DFT); $\Sigma^{(2)}_{\rm phys}={\rm Sgblk}/N_k$ (verified); and the folded route's $Y$ **does not** carry the $1/N_k$, so the ratio target for the §4 V1 acceptance is $1/N_k$, not $1$. **Inside the solver everything is physical** ($S^{\rm phys}$, $\Delta V_{\rm phys}$, $\Sigma_{\rm phys}=-\langle S^{\rm phys}_a|X_b\rangle$, with no stray factors); the $N_k$ is multiplied back only when writing the `.dat`, for compatibility with existing post-processing (`H_eff = diag(ε) + (M+Σ)/N_k`).

## 3. Piece B v2: a cube-anchored `apply_dV`

**Preparation (one-off, ~30 lines):** rearrange `V_colin` (240×240×300) into `dV_blk(nnr_p, 36)` (a 40×40×300 unit-cell grid × 36 cells; startup asserts on commensurability, $240=6\times40$, $300=300$, k-grid $12\times12\times1$); plus a phase table $E(k,m)=e^{i\mathbf k\cdot\mathbf R_m}$ ($144\times36$). The anchor is aligned by construction: cube origin = unit-cell origin, $m=0..5$, and the composed intra- plus inter-cell phase equals `build_V_folded`'s `irx*d1` convention point by point.

**Each application (single source vector, ~6–8 s):**

```text
1. per channel:  ũ_k(r) = e^{i k·r} ∘ invfft(X(:,k))            # 144 unit-cell FFTs + intra-cell phase
2. expand:       Psi_blk(nnr,36) = ũ(nnr,144) · E(144,36)        # ZGEMM -> the true wavefunction on the cube region
3. multiply:     W_blk = dV_blk ∘ Psi_blk                        # pointwise (17M)
4. fold back:    ṽ(nnr,144) = (1/N_k) · W_blk · E^H              # ZGEMM
5. per channel:  Y(:,k') = gather_G[ fwfft( e^{-i k'·r} ∘ ṽ_k' ) ] # 144 FFTs, igk_all mapping
6. nonlocal:     Y(k') += Σ_I |β_I(k')> D_I · (1/N_k) Σ_k <β_I(k)|X(k)>   # KB is separable; reuse get_betavkb/make_coeff
7. apply_Qproj per channel
```

**The nonlocal term must be there** (it is part of $\Delta V$, and the source in (7a) already contains nl; v1's folded comparison against the local part only was an acceptance convention, not a physical one).

## 4. The acceptance ladder: V0 (element-by-element against the $M$ block) is the foundation

| level | test | criterion | cost |
|---|---|---|---|
| **V0 (foundation)** | `apply_dV` (unit source $=|b\rangle$ in its home-$k$ channel) contracted against all $\langle a|$ | $=M_{\rm blk}(a,b)/N_k$ **element by element** (sweeping a few $b$ across whole columns verifies local + nonlocal, phases, anchoring and normalization all at once) | seconds, in-run against the $M_{\rm blk}$ already in memory |
| V1 | against the folded $\Sigma^{(3)}$ (nl off, local against local) | ratio $=1/N_k$ **exactly** | reuses the `do_sigma3` harness |
| V2 | solver regression: switch off the $\Delta V$ term in the matvec | $\Sigma$ column by column $={\rm Sgblk}/N_k$ (direct regression against the stored 2nd-order block) | small subset |
| V3 | physical acceptance: $H_{\rm eff}$ after the full solve | **$e$ returns to $+1.2$–$1.4$**, $a_1$ hugs the VBM; compared against explicit / Koster–Slater | full run |
| V4 | $\max\lvert\Sigma-\Sigma^\dagger\rvert$; spot-check two values of $\omega_0$ | Hermitian; sensible $\omega_0$ behaviour | free |

If V0 does not pass, nothing downstream moves — phase/anchor bugs are laid bare element by element at V0, instead of showing up as a single wrong scalar the way they did in v1.

## 5. Solver and parallel layout (two substantive revisions)

**(1) Source-parallel, with a zero-communication matvec.** The 2nd-order code is pool-per-channel; the all-order matvec couples all channels, so keeping channel parallelism would require an `mp_sum` of a ~280 MB cube field every step — unacceptable. Instead: **each rank holds the complete all-channel vectors for $N_A/36\approx44$ sources** ($\approx15$ MB per source), with `h_psi` (per channel, batched over sources) and `apply_dV` both entirely rank-local, and a single `mp_sum` at the end to assemble $\Sigma$. Sources are naturally parallel ⇒ linear speedup across nodes. New plumbing: `hpsi_setup_globalk(kg)` — overwrite the local `igk_k` slot with the already-gathered `igk_all/xkc`, compute `g2kin` by hand, and fill `vkb` via `get_betavkb` (with save/restore, ~40 lines), letting any rank run `h_psi` for any global $k$. Startup self-check: every rank runs the $\langle\psi|H_0|\psi\rangle=\varepsilon$ gate once for **all 144 channels** (the global version of the existing gate), which immediately exposes any igk/vkb/g2kin mismatch.

**(2) Operator indefiniteness: no more betting on positive definiteness — measure the spectrum first.** The earlier assertion that "rest lies above $\omega_0$ ⇒ $A\succ0$ ⇒ real CG" is **not rigorous**, with two risks:

- **The deep bands 1–6** (Mo 4s4p semicore at $\sim-30$ to $-60$ eV, S 3s at $\sim-12.5$ eV) lie below $\omega_0$, so $(H_0-\omega_0)$ is **negative** on them — meaning the existing 2nd-order operator is already mildly indefinite. It has always converged and passed every check, which is **survival in practice, not a theorem**.
- Adding $\Delta V$ (the vacancy core region is $+11.8$ Ry) can push eigenvalues of $Q(H_0{+}\Delta V)Q$ close to $\omega_0$ (rest-space resonances) ⇒ near-singular.

Strategy: keep the verified **deflated block-CG skeleton plus breakdown monitoring** (report the moment $p^\dagger Ap$ changes sign); have the smoke run **measure** the extreme Ritz values of $A|_Q$ from the Lanczos coefficients along the way; if negative directions cause trouble ⇒ switch to **MINRES** (the proper answer for symmetric indefinite, ~100 lines, same storage); if near-resonant ill-conditioning appears ⇒ a complex solve at $\omega_0+i\eta$. The preconditioner stays per-channel Jacobi $1/\max(g_{\rm kin}^2,1)$. $\alpha=2(\omega_0-{\rm win_{min}})$ is unchanged.

## 6. Cost and memory (anchored to measurements)

| item | anchor | value |
|---|---|---|
| `apply_dV` per source·iteration | 288 unit-cell FFTs $\approx2$ s + 2 ZGEMMs ($4.8\times10^5\times144\times36$) $\approx4$ s + nl $\approx1$ s | **~6–8 s** |
| `h_psi` per source·iteration (batch-amortized) | the per-channel cost of the 2nd-order block (1 h 59 m) | ~1 s |
| $n_{\rm iter}$ | to be measured in the smoke run ($\Delta V$ core is $+11.8$ Ry and the preconditioner only covers kinetic energy, so possibly higher than 2nd-order) | ~100 (assumed) |
| **full block, 1584 sources** | $1584\times100\times7\,{\rm s}/36\ {\rm ranks}$ | **$\approx$8.6 h (1 node)**; zero communication ⇒ $\div N$ (4 nodes $\approx$2 h) |

Memory per rank: the $X$ batch (8 sources × 5 CG vectors × 15 MB) $\approx0.6$ GB + cube buffers $2\times0.28$ GB + `evc_act_all` 165 MB $\approx$ **1.5 GB** ✓. Merging ZGEMMs over small source batches and caching `vkb` per channel per iteration can save another 1.5–2×.

## 7. Phase plan, risks, and where this sits

| phase | content | acceptance | effort |
|---|---|---|---|
| **P-I** | `dV_blk` rearrangement + `apply_dV` v2 (with nl) + the V0/V1 harness | V0 passes element by element; V1 ratio $=1/N_k$ | ~1–2 days |
| P-II | `hpsi_setup_globalk` + matvec + the all-channel gate | V2 regression $={\rm Sgblk}/N_k$ | ~1 day |
| P-III | block solve (deflation + monitoring) → smoke (288 columns, bands 13–14) → all 1584 columns @ $\omega_0{=}$VBM, 2–4 nodes | Hermiticity; measured $n_{\rm iter}$ / spectrum | ~1 day + 2–4 h of compute |
| P-IV | $H_{\rm eff}$ levels + write-up | $e$ back to $+1.2$–$1.4$, compared against explicit/KS | half a day |

**Risk table:** ① phases/anchoring (v1's killer) → nailed down element by element at V0; ② indefiniteness / iteration count → measure the spectrum first, with MINRES / $+i\eta$ in reserve; ③ a missing or double-counted nl term → V0 includes nl; ④ module-state side effects from running `h_psi` at a global $k$ → save/restore plus the all-channel gate; ⑤ cube buffer memory → bounded as long as source batches are $\le$4.

**Where this sits.** If all you want is defect levels, [Koster–Slater](koster-slater.html) (minutes) and explicit (46 min) already give the same $a_1{+}e$ — the value of this plan is that it **treats the root cause** of the rest dressing (all-order $\Delta V_{QQ}$, removing second order's over-screening and divergence) and produces a spectral-function-grade production object: one static $\Sigma(\omega_0)$ block costs ~half a day per node, and when the full $\Sigma(\omega)$ is wanted later, P-III simply repeats per $\omega$ (amortizing linearly across nodes).

## 8. Execution checklist: task breakdown and acceptance criteria

§7's phase table expanded into a task breakdown with acceptance criteria (numbers #11–#26 correspond to the work tracker; an item is ticked when its *acceptance criterion* is met, **not** when the code is written). Dependency chain: $11\to12\to13\to14$; $15$ can run in parallel with P-I; $\{14,15\}\to16\to17\to18\to19\to20\to21\to\{22,23\}$; $22\to24$; #25/#26 are an independent side branch. **The critical path is $11\to12\to13$ (V0)** — once V0 passes, the rest is wiring and machine time.

### P-I Piece B (`apply_dV`) + acceptance — the only new physics machinery (~1–2 days)

| # | task | acceptance criterion |
|---|---|---|
| 11 | `dV_blk` preparation: rearrange `V_colin` into $(40{\times}40{\times}300,\,36\ \text{cells})$ blocks + phase table $E(k,m)=e^{i\mathbf k\cdot\mathbf R_m}$ + commensurability asserts ($240=6\times40$ etc.) | any single cell agrees point by point with `V_colin` in its original order |
| 12 | `apply_dV` v2: intra-cell phase $\to$ expansion ZGEMM $\to$ multiply by `dV_blk` $\to$ fold back (with the $1/N_k$) $\to$ per-channel FFT/gather $\to$ nonlocal (with an nl-off switch) $\to$ Qproj | compiles; acceptance via #13/#14 |
| 13 | **V0, the foundation**: a unit source ($=|b\rangle$ in its home-$k$ channel) through `apply_dV`, contracted against all $\langle a|$, for ~8 values of $b$ spanning bands and $k$ | **element by element** $=M_{\rm blk}(a,b)/N_k$, relative error ~$10^{-10}$; on failure go back to #12 and freeze everything downstream |
| 14 | V1: nl-off, applied to the retained $\chi^0$, contracted against the folded $\Sigma^{(3)}$ | ratio $=1/N_k=6.944\times10^{-3}$ **exactly** (v1's failing value was $4.4\times10^{-5}$) |

### P-II Global matvec (~1 day)

| # | task | acceptance criterion |
|---|---|---|
| 15 | `hpsi_setup_globalk`: overwrite the igk slot + compute $g^2_{\rm kin}$ by hand + fill vkb via `get_betavkb`, with save/restore (**can start in parallel with P-I**) | the $\langle\psi_a|H_0|\psi_a\rangle=\varepsilon_a$ gate at ~$10^{-9}$ eV on every rank for all 144 channels |
| 16 | assemble the all-order matvec: $A=Q(H_0{+}\Delta V{-}\omega_0)Q+\alpha P$, source batches $\le8$/rank, cube buffers reused | $\langle Y|AX\rangle=\langle AY|X\rangle^*$ to ~$10^{-12}$ on random Q-vectors |
| 17 | V2 regression: switch off the $\Delta V$ term in the matvec and solve $\Sigma$ on a small source subset | column by column $={\rm Sgblk}/N_k$ to the CG threshold — the end-to-end regression of the new layout (source-parallel + global $k$) |

### P-III Solver + runs (~1 day + machine time)

| # | task | acceptance criterion |
|---|---|---|
| 18 | `solve_feshbach_block`: block-CG + deflation + **breakdown monitoring** (report the moment $p^\dagger Ap$ changes sign) + Ritz spectrum diagnostics; MINRES / $\omega_0{+}i\eta$ only if the monitor fires | convergence + a spectrum report (definite or not) |
| 19 | `do_feshbach` driver: source construction, source-parallel distribution ($\approx44$ sources/rank), assembly of $\Sigma_{ab}=-\langle S^{\rm phys}_a|X_b\rangle$, writing the `.dat` (multiplying $N_k$ back at the boundary, same format) | memory of ~1.5 GB/rank confirmed |
| 20 | smoke run: 288 columns (bands 13–14, all $k$), 1 node | measured $n_{\rm iter}$ (assumed ~100; $>300$ means changing the preconditioner or going to $+i\eta$), the Ritz spectrum decides the solver, wall clock calibrated |
| 21 | production run: all 1584 columns @ $\omega_0{=}$VBM, 2–4 **exclusive** nodes (36 ranks/node, unshared) | $\max\lvert\Sigma-\Sigma^\dagger\rvert$ before Hermitization at the same order as 2nd-order's $5.9\times10^{-14}$; ~2–4 h |

### P-IV Physical acceptance + write-up (~1 day)

| # | task | acceptance criterion |
|---|---|---|
| 22 | $H_{\rm eff}$ levels + a fourth ΔDOS panel (same mesh and $\eta$ as Fig 15) | **$e$ returns from second order's $+0.36$ to $+1.2$–$1.4$** (explicit $+1.35$ / DFT $+1.19$), $a_1$ hugs the VBM, degeneracies $1{+}2$ — the final acceptance for the whole line |
| 23 | (optional) a second reference production run at $\omega_0'\approx-4.6$ + assessment of the static-$\omega_0$ error | the $e$ values at the two $\omega_0$ differ by $<0.1$ eV $\Rightarrow$ the static approximation is settled |
| 24 | write-up: result figures, complete the Fig 13/15 narrative with the all-order conclusion, turn this page's Drafted numbers into measured ones, flip the catalog to prod, update memory, push | build self-checks pass and the published set is clean |

### Side branch: Koster–Slater (independent, can run at any time in parallel)

| # | task | acceptance criterion |
|---|---|---|
| 25 | host GF band convergence: enlarge the Wannier band space and run the $e(N_{\rm bands})$ convergence curve | $e$ moves from $+1.50\to+1.35$/$+1.19$; Krein–Friedel count $2.06\to3$ |
| 26 | (optional) $C_{3v}$ characters: use the $\sigma_v$ character to fix $a_1$ vs $a_2$ rigorously + irrep factorization | removes the "by analogy" caveat on the [KS page](koster-slater.html); incidentally a $6$–$12\times$ symmetry speedup |

**Status (2026-06-09):** all pending; currently unblocked and startable immediately are **#11 and #15** (the parallel entry points of the Feshbach line) and **#25, #26** (the side branch).
