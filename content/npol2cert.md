# npol=2 certification: the verdict on the doubling test, and its provenance

## 0. Summary

The npol=2 (spinor) port of the MODE C downfolding chain failed to reproduce "levels exactly ×2" in the doubling
test (noncolin=T, lspinorb=F, scalar pseudopotentials), which triggered a full forensic investigation. Two acts,
two conclusions. Act one: **EDT's npol=2 implementation is innocent** — all three independent gauge-invariant
criteria pass at machine precision. Act two (a reversal, after pair-by-pair instrumented dumps): **EDI's compute
engine is innocent too (exact pair by pair to $10^{-12}$); the real culprit is the file interface** — EDI-direct's
.bin writes $(\mathrm{bra},\mathrm{ket})$ within each $(k_i,k_f)$ block, while EDT's read convention is
$(\mathrm{ket},\mathrm{bra})$ — **a pure band-axis transpose**, with the numerical verdict
$\max|E_\mathrm{EDI}-T_\mathrm{EDT}^{T(\mathrm{bands})}|=8.7\times10^{-14}$.
Fix: transpose on the EDI write side (format version 20260702→20260809); on the EDT side, the born block $M_{AA}$
was ported into EDT proper (npol=2), making the SOC production line self-sufficient. Final verdict: spinor-pair
degeneracy **0.0000 meV**, cross-host VBM-gated agreement ≤4.9 meV (the host-twin noise floor), chain-vs-born
closure $4.6\times10^{-11}$, and after the fix EDI-vs-EDT on the same wavefunctions agrees over the full matrix
at $\sim10^{-13}$.

## 1. The doubling test and the initial failure

Construction of the doubling test: the same MoS$_2$ V$_\mathrm{S}$ 6×6 system, with the host recomputed using
noncolin=T, lspinorb=F and **scalar pseudopotentials** (d66nc). Then $H$ is spin-diagonal with no magnetization,
so the physics is identical to the scalar calculation: every level of the downfolded effective Hamiltonian should
appear exactly twice. The first run (with EDI-direct supplying the born block $M_{AA}$) FAILED: the levels neither
paired up nor corresponded, and the chain's operator identity test $\max|P^\dagger(\Delta V\,\Psi)-M|$ reached
$8.7\times10^{-1}$ (the scalar case is of order $10^{-14}$).

## 2. Three-way decomposition: try the local arm and the KB arm separately

Two probe hooks split apply_dV apart: `EDT_NOKB=1` disables the KB (nonlocal) channel to leave the pure local
projection, and an independently written numpy plane-wave integrator (the referee — a zero-dependency
reimplementation: zero-padded ifftn, $e^{iq\cdot r}$ phases, sum over both spinor components) reconciles all 160
complex matrix elements one by one:

| audited object | criterion | result |
|---|---|---|
| EDT spinor local fold | $\|P_\mathrm{nokb}\|/\|\mathrm{referee}\|$ | 52/52 element ratios **1.000** (all within 3%) |
| EDT spinor KB | $\|P_\mathrm{full}-P_\mathrm{nokb}\|/\|M_\mathrm{EDI}-\mathrm{referee}\|$ | median 0.36, p10=0, spread 0.13–1.3 |

The local arm was acquitted on the spot (and the referee and EDT cross-verify each other). The KB arm's apparent
"evidence of guilt" was later shown to be a wrongful conviction — the denominator
$M_\mathrm{EDI}-\mathrm{referee}$ was itself contaminated.

## 3. The key evidence: mixed spinors and degenerate-pair parallelism

d66nc's eigenstates are **heavily mixed spinors** (purity 0.51–0.99: Davidson mixes arbitrarily in SU(2) within a
degenerate pair, so a pure-spinor assumption does not hold). For a spin-diagonal operator with no SO, the true
matrix elements satisfy

$$\langle m|\mathcal{O}|n\rangle = \langle\chi_m|\chi_n\rangle\,\mathcal{O}_{\mathrm{orb}}(m,n),$$

so for a degenerate pair $(m_1,m_2)$ (same orbital, orthogonal spinors) the **KB 2-vector
$[K(m_1,n),K(m_2,n)]$ must be complex-proportional (parallel) to the local 2-vector
$[L(m_1,n),L(m_2,n)]$** (both $\propto[\langle\chi_{m_1}|\chi_n\rangle,\langle\chi_{m_2}|\chi_n\rangle]$).
This is the ultimate gauge-invariant criterion — completely immune to per-k unitaries and to degenerate mixing:

| KB under audit | parallelism $|\cos|$ | pass rate (>0.99) |
|---|---|---|
| EDT ($P_\mathrm{full}-P_\mathrm{nokb}$) | median **1.0000**, p10 1.000 | **25/25** |
| EDI ($M_\mathrm{EDI}-\mathrm{referee}$) | median 0.977, p10 0.759 | 44% |

At the same time the "EDI total residual" $M_\mathrm{EDI}^{\rm file}-(\mathrm{referee}+K_\mathrm{EDT})$ reached 0.50
(against $\max|M|=0.42$) and correlated with the local amplitude (0.69) — at that point the chain of evidence
pointed at "EDI matrix elements, local and KB alike, are untrustworthy", and on that basis the born block was
ported into EDT (§4). **But the mechanism behind that conviction was later overturned** (§7): pair-by-pair
instrumented dumps proved the EDI engine exact pair by pair, and the fault lay in the file interface's band-axis
convention. (An audit incidentally turned up a $\sigma_1$-only contraction in the `ed_coarse_full_q` noncolin
branch, a genuine code bug now fixed, but it lives on the interpolation path and was not in this case's
propagation chain.)

**Methodological lessons (on the record)**: (1) the earlier singular-value "exoneration" of EDI (median 6.5e-3)
held yet was untrustworthy for one root reason — **singular values are invariant under transpose**, so SV-class
criteria are natively blind to this class of interface error; (2) full-matrix Hermiticity makes the amplitude
spectra of $|M^T|$ and $|M|$ identical, and "median ratio 1.0000 with a large p10/p90 spread" is precisely the
fingerprint of a transpose; (3) **always look at the whole distribution, and always add an element-by-element
complex-level reconciliation.**

## 4. The fix: born block into EDT, the whole noncolin path de-EDI-ed

`vtilde_block_mpi` (the born path, `born_only=.true.`) was ported to npol=2: spin dimensions added to the source
vector $S$ and to the bec/coeff arrays, a component loop added in the fold ($\Delta V_\mathrm{loc}$ is
spin-diagonal), KB routed through the same `make_coeff_nc` as apply_dV (dvan spin-diagonal / dvan_so 2×2), and the
bra side contracted component by component; npol=1 is bit-for-bit unchanged. The verification ladder (kbval):

| gate | result |
|---|---|
| scalar born regression (r9 vs archive) | **bit-identical** (8,296,152 B) |
| nc born (40 bands × 36 k, $N_A=1440$) | Hermiticity $1.1\times10^{-14}$ |
| chain-vs-born closure unit test | $\mathbf{4.6\times10^{-11}}$ (0.87 against the EDI edmat) |
| SV doubling audit, **all** 1260 off-diagonal blocks | median 3.5e-3, p90 7.8e-3, max 1.1e-2 (the level of a two-SCF host difference, no bad tail) |

## 5. Final verdict: a two-gate framework

The original "exactly ×2 to within 0.2 meV" implicitly assumed **the same Hamiltonian**. In reality the scalar host
(dout66f) and the noncolin host (d66nc) are two independent SCF runs, and QE's noncolin XC (the $m\to0$ path) makes
the raw eigenvalues differ by 3–52 meV and the VBM by 10.8 meV — 0.2 meV is unreachable in principle across
different hosts. The correct framework:

| gate | nature | result |
|---|---|---|
| spinor-pair degeneracy $\max|e_{2j}-e_{2j-1}|$ | same host, **exact** | **0.0000 meV — PASS** |
| scalar-vs-nc after gating on each VBM | cross-host, limited by twin noise | ≤4.86 meV — PASS |

Together with the machine-precision internal evidence of §2–§4, npol=2 EDT (chain + born) reaches production
certification. Acceptance and forensic scripts live in the qe-edt repository under `post/`: `doubling_accept.py`,
`parallel_test.py`, `analyze_split.py`, `svdouble2.py`.

## 6. Direct consequences for SOC production

- SOC (lspinorb=T) chains always use **EDT's own born edmat** (`do_full_block=.true., born_only=.true.`);
  feeding them an EDI-direct noncolin/SOC edmat is forbidden.
- An incidental find: the old SOC O$_\mathrm{S}$ supercell used **unrelaxed geometry** (the O substitutional $z$
  off by 1.0 bohr), inconsistent with the non-SOC production line (relaxed) — all four SOC supercells (clean,
  O$_\mathrm{S}$, V$_\mathrm{S}$, Se$_\mathrm{S}$) are queued for recomputation on the campaign geometry.
- 50 vs 100 Ry convergence check: the SOC splittings themselves are converged at 50 Ry (K-VB 149.0 meV, K-CB
  3.0 meV, differing by <0.1 meV), but absolute energies in the active window drift by up to 16 meV → the SOC stack
  stays at 100 Ry.

## 7. Act two: pair-by-pair instrumented dumps and the interface culprit

Fitting EDI with an `EDI_DBG_PAIR` probe (dumping the pre-assembly $m_\mathrm{loc}$ and $m_\mathrm{nl}$ matrices
for a specified $(k_i,k_f)$ pair) allowed a three-way reconciliation against the referee and the already-certified
EDT KB:

| reconciliation (pair $(k_1,k_{29})$, 80 elements) | result |
|---|---|
| EDI $m_\mathrm{loc}$ vs independent referee | $6.9\times10^{-12}$ |
| EDI $m_\mathrm{nl}$ vs EDT KB (certified) | $9.5\times10^{-14}$ |
| but assembled .bin vs the EDT-born full matrix | max 1.68 |

Every number right, every placement wrong → the fault is in assembly/writing. A single sweep over 7 axis
hypotheses:

$$\max\bigl|E_\mathrm{EDI} - T_\mathrm{EDT}^{T(\mathrm{bands})}\bigr| = 8.7\times10^{-14}$$

— a pure band-axis transpose within the block (no $k$ exchange, no conjugation). Fixed by transposing on the EDI
write side, with format version 20260809 (edi_v8d); **every 20260702-version EDI-direct .bin carries the band
transpose** and must be regenerated or transposed on read (files with that same version number produced by conv.py
follow the EDT convention and are unaffected). Fixed in the same batch: B8 ($q_\mathrm{cryst}$ now taken in the
NSCF representation, preventing spurious $e^{i\Delta G\cdot r}$ phases) and the full_q $\sigma_1$-only bug. The
final closure after the fix (edi_v8d, same set of wavefunctions):

| closure | before the fix | after |
|---|---|---|
| noncolin: EDI .bin vs EDT-born full matrix (40 bands × 36 k) | max 1.68 | **max $8.7\times10^{-14}$, median $2.1\times10^{-16}$** |
| scalar: EDT chain unit test vs EDI edmat (cross-code) | 1.677 | $\mathbf{1.86\times10^{-14}}$ |

Two independent implementations verifying each other to machine precision — the strongest cross-code closure of
the whole campaign.

**Overall lesson**: when two codes disagree, the first move is not to audit the engines but to **separate "engine"
from "interface" with pair-by-pair instrumented dumps** — two individually correct engines can perfectly well
convict each other across a single transpose.
