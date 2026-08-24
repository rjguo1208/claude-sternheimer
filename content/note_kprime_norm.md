# Implementation Note — the $k'$-sum normalization convention

*Why the rest-space sum in the Layer-1 downfolding carries a $1/N_k$ factor, and how that settles the "channel-by-channel sum is far too large" problem in one stroke. EDT P2 → P3; see the [Implementation Plan](plan.html).*

> **Conclusion (settled):** the discrete sum over $k'$ *is* the discretization of a BZ integral, so it carries a **$1/N_k$** (BZ measure) factor. That turns the naive sum's $\sim-70$ Ry into a physical $\sim-0.5$ Ry. Furthermore, as long as $\Delta V$ is converged inside the supercell, the matrix element $M$ (and hence $\tilde V$) is **independent of supercell size**, so there is **no** extra $N_{\rm sc}$ factor beyond $1/N_k$ (§4). Below: the derivation, its reconciliation with "channel-by-channel closure", the supercell independence, and the numerical checks used to confirm it.

## 0. The quantity in question

The rest self-energy (the second-order term of the Feshbach downfolding):

$$
\Sigma_{mn}(\omega_0)=\big\langle m k_i\big|\,\Delta V\,G^R(\omega_0)\,\Delta V\,\big|n k_i\big\rangle
= \frac{1}{N_k}\sum_{k'}\sum_{r\in R(k')}\frac{M_{m,\,rk'}\;M_{rk',\,n}}{\omega_0-\varepsilon_{rk'}},
\qquad M_{m,rk'}\equiv\langle m k_i|\Delta V|r k'\rangle .
$$

That $1/N_k$ is the measure factor argued for below ($N_k$ = number of rest BZ mesh points). The downfolded potential $\tilde V_{mn}=M_{mn}+\Sigma_{mn}$ has the same dimension and normalization as $M$, and the Born limit reduces to $\tilde V\to M$.

## 1. Why $1/N_k$: the $k'$ sum *is* a BZ integral

In an infinite crystal the resolution of the identity is a BZ integral:

$$
\mathbb 1=\sum_n\frac{\Omega_{\rm cell}}{(2\pi)^3}\int_{BZ}\!d^3k\;|\psi_{nk}\rangle\langle\psi_{nk}| .
$$

Discretizing on an $N_k$-point Monkhorst–Pack mesh, $\dfrac{\Omega_{\rm cell}}{(2\pi)^3}\displaystyle\int_{BZ}d^3k\to\dfrac1{N_k}\sum_k$ (because $V_{BZ}=(2\pi)^3/\Omega_{\rm cell}$), and with QE's **per-unit-cell normalized** states ($\langle\psi_{nk}|\psi_{nk}\rangle_{\rm cell}=1$) this gives

$$
\boxed{\;\mathbb 1=\frac{1}{N_k}\sum_{nk}|\psi_{nk}\rangle\langle\psi_{nk}|\;}
$$

so both the resolvent and the rest Green function carry $1/N_k$:

$$
G_0(\omega)=\frac{1}{N_k}\sum_{nk}\frac{|\psi_{nk}\rangle\langle\psi_{nk}|}{\omega-\varepsilon_{nk}},
\qquad
G^R(\omega_0)=\frac{1}{N_k}\sum_{k'}\sum_{r\in R}\frac{|rk'\rangle\langle rk'|}{\omega_0-\varepsilon_{rk'}} .
$$

**Consistency**: Layer-2's $G^A(\omega)$ carries $1/N_k$ in the same way; the convergence of $\frac1{N_k}\sum_{k'}\to\int_{BZ}$ under rest-mesh refinement is the quantitative version of the earlier correction that "rest must cover the whole BZ"; and the $\frac1{N_k}\sum_{k'}$ in the EDI golden rule (the final-state BZ average) is a **separate but identically-sourced** $1/N_k$. In short: **every internal $k$ sum carries a $1/N_k$.**

## 2. Reconciling with "channel-by-channel closure" (two different things)

At fixed $k'$, summing over **all bands** is Parseval within that $k'$ sector, and carries **no factor**:

$$
\sum_{n'\,(\text{all})}\big|\langle n'k'|s\rangle\big|^2=\langle s(k')|s(k')\rangle_{\rm cell}=\sum_G|s(k',G)|^2,
\qquad |s\rangle=\Delta V|n k_i\rangle .
$$

Summing that (factor-free) per-channel result over $k'$ **with the $1/N_k$** gives a finite, unit-cell-scale quantity:

$$
\frac{1}{N_k}\sum_{k'}\sum_{n'\,(\text{all})}|M_{n,n'k'}|^2=\frac{1}{N_k}\sum_{k'}\langle s(k')|s(k')\rangle=\big\langle n k_i\big|\Delta V^{2}\big|n k_i\big\rangle .
$$

So **the two statements do not conflict**: the band sum within one channel is sector Parseval (no factor), while the sum across $k'$ is a BZ measure (with $1/N_k$). Conflating them is what led to the earlier, incorrect conclusion of "no $1/N_k$".

## 3. Orders of magnitude: $-70$ Ry → $-0.5$ Ry

P2's Sternheimer solve gives a **per-channel** $\Delta\tilde V_{\rm chan}(k')\sim0.3$–$0.7$ Ry. Naively adding 144 channels gives $\sim-70$ Ry (unphysical). With the $1/N_k$:

$$
\Sigma_{nn}=\frac{1}{N_k}\sum_{k'}\Delta\tilde V_{\rm chan}(k')\approx\frac{-70}{144}\approx-0.5\ \text{Ry},
$$

the same order as the Born self-energy $M_{nn}\sim0.7$ Ry — a **large but physical** correction (a vacancy is a strong scatterer, Born *ought* to be substantially corrected, and this is exactly where the $T$-matrix earns its keep).

## 4. Supercell independence ⇒ no $N_{\rm sc}$ factor

**Physical requirement:** as long as $\Delta V$ is **converged** inside the supercell (localized, decayed to zero before the boundary), the electron–defect matrix element $M=\langle\psi_{nk_i}|\Delta V|\psi_{mk_f}\rangle$ is **independent of supercell size** — a larger supercell only adds bulk/vacuum region where $\Delta V=0$, which cannot change an integral that only receives contributions near the defect.

**The implementation agrees:** EDI's $M=\frac1{N_{\rm nnr}}\sum_{r}u^*\,V^q_{\rm folded}\,u$ uses the **unit-cell** measure $\frac1{N_{\rm nnr}}$ (not the supercell measure $\frac1{N_{\rm nnr}^{\rm SC}}$), acting on the $V^q_{\rm folded}$ that folds the localized $\Delta V$ back into the unit cell. Locality ⇒ the sum only receives contributions in the defect region (fixed) ⇒ $M$ **does not scale with $N_{\rm sc}$** (T1 already showed $M$ agrees with EDI to $2.3\times10^{-13}$). Single-defect concentration is handled by $n_d$ in the golden rule and **does not enter $M$**.

**Corollary:** $\Sigma=\frac1{N_k}\sum_{k'}\sum_r MM/(\omega_0-\varepsilon)$ is built from supercell-independent $M$, so $\Sigma$ and $\tilde V$ are **supercell-independent** too, and the $\Sigma_{nn}\approx-0.5$ Ry estimated above is a **per-defect physical quantity**. **Hence there is no extra $N_{\rm sc}$ factor beyond $1/N_k$.**

**Two independent convergence axes (do not conflate them):** (i) **supercell size** — decides whether $\Delta V$ (and hence $M$) is converged and supercell-independent; (ii) **rest BZ mesh $N_k$** — the integration convergence of $\frac1{N_k}\sum_{k'}\to\int_{BZ}$ ($N_k$ may be denser than the supercell folding allows, by taking $\Delta V(q)$ from a real-space folding at arbitrary $q$).

## 5. Confirmation checks (the P3 starting point)

$1/N_k$ is fixed by the BZ measure and $M$ is supercell-independent (§4), so the normalization is fully determined. The following checks **confirm the implementation is correct**:

- **Per-channel closure**: at a single $k'$, $\sum_{n'=1}^{150}|M_{n,n'k'}|^2\overset?=\sum_G|s(k',G)|^2$ (the shortfall = high-band incompleteness, and should be consistent with the Sternheimer high-band tail).
- **Total sum $=\langle\Delta V^2\rangle$**: $\dfrac1{N_k}\sum_{k'}\sum_{n'}|M|^2\overset?=\dfrac1\Omega\int_{\rm cell}|\psi_{nk_i}|^2\,\Delta V^2$; both sides should be supercell-independent.
- **Supercell independence**: given two supercell sizes, $M$ (and $\Sigma$) should be unchanged; EDI transport already verifies this implicitly (converged mobility is supercell-independent).
- **Born-limit anchor (the gold standard)**: switch off the rest/active resummation ⇒ $\tilde V\to M$, $T_{PP}\to M$, run transport, and it **must reproduce the EDI mobility bit for bit** (including $n_d$ and $\frac1{N_k}$).

> **Summary**: the normalization is fully settled — the $k'$ sum carries the BZ measure $1/N_k$; $M$ is **supercell-independent** because $\Delta V$ is localized (no $N_{\rm sc}$ factor); the per-channel band sum carries no factor, and $\frac1{N_k}\sum_{k'}\sum_{n'}|M|^2=\langle\Delta V^2\rangle$. Two convergence axes: supercell size ($\Delta V$) and rest BZ mesh ($N_k$). The overall anchor is the Born limit reproducing the EDI mobility.
