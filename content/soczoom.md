# SOC defect spectral functions: K-valley band-edge zooms (first batch)

## 0. Summary

Electron–defect spectral functions $A(k,\omega)$ **including spin–orbit coupling**: monolayer MoS$_2$, three defects O$_\mathrm{S}$ /
Se$_\mathrm{S}$ / V$_\mathrm{S}$ (relaxed geometry), a 12×12 coarse mesh × 22-spinor-band chain
(bands 13–34) plus a **22-orbital Wannier vertex frame** (the same wannierization, no disentanglement),
$n_d = 10^{12}\,\mathrm{cm^{-2}}$, $\eta=5$ meV, the KZOOM path ($K-0.3\,\overline{\Gamma K}
\to K+0.45\,\overline{KM}$, 721 points), $N_f=1200$. The CBM window shows the **conduction-band SOC doublet**
(3.0 meV splitting at K, opening along $K\to M$); the VBM window shows the **149 meV K-point valence spin splitting**
with the defect flat-band multiplet between the two. In the full-path window, V$_\mathrm{S}$'s deep in-gap S-p doublet (+1.12/+1.17)
appears right across the path. Spectral functions come from the fast route (vertex factorization + host $G_0$ spectral binning, `kpath_fast.py`):
all 9 windows in ~30 minutes, with in-gate verification that the factorization is exact (1e-15) and the binning is 3e-5 mid-gap / 7e-4 at the band edge.

## 1. Figures

**Full-path overview ($\Gamma$–M–K–$\Gamma$, full energy window; 2026-08-12, 22-orbital Wannier frame
[bands 13–34, no disentanglement]). In the V$_\mathrm{S}$ panel the **deep in-gap S-p doublet at +1.12/+1.17 eV**
shows up across the whole path (consistent with the dilute-DOS levels of §4b; $\max|T|$ in $T(\omega)$ is enhanced
×200 at both levels), while the O$_\mathrm{S}$/Se$_\mathrm{S}$ gaps are clean — the three-defect comparison is
itself the physical criterion. Colour-scale floor 0.03: an isolated defect flat band at $n_d=8.8\times10^{-4}$/cell
simply carries little spectral weight ($\max_k A \approx 0.2$), and this is the honest brightness of the dilute limit:**

![SOC22 full path](../assets/kfull_soc22.png)

> Process note: the earlier 10-orbital version was missing this doublet for **two** independent reasons — the
> Mo-d subspace cannot represent it, and the chk→u.mat converter wrote an effectively transposed $U$
> (a unitarity check is blind to that). The latter was diagnosed and fixed with the wannierization identity
> $U^\dagger\,\mathrm{diag}(\varepsilon_k)\,U = H_W(k)$
> (1.5e-5 for the correct orientation, O(4 eV) for the wrong one; the tool now lives in qe-edt as `post/chk2umat.py`).

**Three-defect 22-band K-valley zooms (144k chain + 22-orbital frame, same wannierization):**

![SOC22 CBM zoom](../assets/kzoom_soc22_cbm.png)

![SOC22 VBM zoom](../assets/kzoom_soc22_vbm.png)


**Self-energy colour maps ($\Sigma^{ed}(k,\omega)=n_d\,T(k,k;\omega)$, same KZOOM path and energy windows;
top row Re Tr $\Sigma$ (diverging colour scale), bottom row $-$Im Tr $\Sigma$ (logarithmic)). The continued
fraction is evaluated at $\omega+i\eta$ (retarded $\tilde V$): a finite chain discretizes the dressed rest
continuum into real-axis Ritz poles, and evaluating at real $\omega$ makes those in-window poles appear as
zero-width needles (diagnostic: $\max|S_m|$ = 71/63/161 at the poles, at defect-dependent energies, i.e. a
chain-side artifact rather than a host vHs). With $+i\eta$ the needles become physical Lorentzians of width
~$\eta$, and in-gap levels shift only by $O(\eta^2/\Delta)$:**

![Sigma CBM zoom](../assets/sigzoom_soc22_cbm.png)

![Sigma VBM zoom](../assets/sigzoom_soc22_vbm.png)

Reading the maps: $T(k,k)$ is nearly flat along the path — the signature of a short-ranged point-defect
potential. In the VBM window, V$_\mathrm{S}$'s shallow +0.03 in-gap state is a textbook resonance (a bright
$-$Im line, Re crossing zero and changing sign, dispersive lobes at $\pm$300 meV), while Se$_\mathrm{S}$ shows
almost no feature (isovalent). In the CBM window, O$_\mathrm{S}$'s Re$\,\Sigma<0$ is level attraction from its
+1.545 resonance below, and V$_\mathrm{S}$'s Re$\,\Sigma>0$ growing into the gap is repulsion from the deep
doublet (+1.12/+1.17) below; $-$Im switches on at the band edge in step (scattering phase space).


**The same figures at $n_d = 10^{13}$ cm$^{-2}$** ($\Sigma = n_d T$ is strictly linear and $T$ is reused from
cache — the colour scale simply scales ×10: V$_\mathrm{S}$'s resonance core reaches Re$\,\Sigma$ = ±3.8 eV and
$-$Im exceeds 1 eV. Note that at this concentration the single-point-defect $T$-matrix (independent-scatterer
assumption) is getting strained — the spectral renormalization near the resonance is no longer perturbative, and
defect–defect interference / ATA corrections enter at this level; the $10^{12}$ version is the quantitatively
reliable regime):

![Sigma CBM zoom nd13](../assets/sigzoom_soc22_cbm_nd13.png)

![Sigma VBM zoom nd13](../assets/sigzoom_soc22_vbm_nd13.png)


**Concentration series: $n_d = 3\times10^{13}$ cm$^{-2}$ (2.64%/cell) — the impurity-band regime**.
The spectral function is already strongly reconstructed: in the CBM window O$_\mathrm{S}$ **pulls the CB bottom
down and splits it** (its +1.545 resonance becomes an impurity band), while V$_\mathrm{S}$ **pushes the CB bottom
up** by ~50 meV — directions that are exactly the amplified realization of the Re$\,\Sigma$ signs in the
$10^{12}$ version. In the VBM window, V$_\mathrm{S}$ anticrosses the valence-band top with the defect level and
pushes the edge up ~0.13 eV, O$_\mathrm{S}$ by +0.05, and Se$_\mathrm{S}$ shifts slightly overall while staying
band-like. **Note**: at this concentration the $n_d T$ independent-scatterer approximation is qualitative only
(defect–defect interference / multiple scattering are not included):

![A CBM zoom nd3e13](../assets/kzoom_soc22_nd3e13_cbm.png)

![A VBM zoom nd3e13](../assets/kzoom_soc22_nd3e13_vbm.png)

![Sigma CBM zoom nd3e13](../assets/sigzoom_soc22_cbm_nd3e13.png)

![Sigma VBM zoom nd3e13](../assets/sigzoom_soc22_vbm_nd3e13.png)

The first-batch 10-band version (O$_\mathrm{S}$+V$_\mathrm{S}$), archived:

![SOC CBM zoom](../assets/kzoom_soc_cbm.png)

![SOC VBM zoom](../assets/kzoom_soc_vbm.png)

The green dashed lines in the V$_\mathrm{S}$ panel of the VBM figure are the **scalar (non-SOC) supercell
ground-truth levels**, for a direct comparison of how SOC rearranges and splits the in-gap defect levels.

## 2. Flat-band detector readings (localized defect features, $E-E_\mathrm{VBM}$/eV)

| window | O$_\mathrm{S}$ | V$_\mathrm{S}$ |
|---|---|---|
| CBM | +1.694 | +1.702, **+1.758** (two features, ~56 meV apart) |
| VBM | −0.253, −0.241, −0.217, −0.181, −0.116, −0.060, −0.020, +0.004 | −0.329, −0.297, −0.253, −0.241, −0.213, −0.181, −0.132, −0.060, −0.012, **+0.053, +0.181** (well into the gap) |

The six scalar V$_\mathrm{S}$ ground-truth levels (−0.298…+0.060) occupy the same range as the SOC multiplet — SOC
splits and rearranges each orbital level according to its total-angular-momentum structure, and V$_\mathrm{S}$
gains a new flat band at +0.181 eV on the gap side.

## 3. Computational route (this batch)

- **born $M_{AA}$**: EDI-direct v8d (after the band-axis transpose fix, certified against EDT cross-code at 1e-14),
  all 40 bands over 144×144 k pairs, 5m36s;
- **MODE C chain**: EDT r11 (batched A2A transpose surgery, 6×6 bit-for-bit regression) / Anvil highmem single node,
  16 pools × 3 threads, n_reorth=1 (0.000000 meV verdict), col_chunk=1, ~45 min/step × 16 steps;
  the chain start-up unit test *is* the cross-code closure for every production run (O$_\mathrm{S}$ 2.730e-11, V$_\mathrm{S}$ 2.547e-11);
- **10-spinor-band Wannier**: soc40 host bands 25–34 form a globally isolated group (net gaps of 6.5 meV / 0.56 eV
  to bands 24/35), so no disentanglement is needed; a w90 3.1 write crash was recovered by extracting from the chk
  directly plus restart=plot (u.mat unitarity 1.9e-10);
- **T cache**: Koster–Slater cluster (RCUT=4, 61 cells, dim 610), only ~5 minutes per window at 10 bands
  (the 11-band scalar production took 2–3.5h) — completed with four windows in parallel on 4 Banff nodes;
- host: 100 Ry (a 60 Ry cross-check: splittings agree to ≤0.03 meV, but absolute levels are off by 7.32 meV, over the 5 meV gate).

## 4. Supercell ground-truth gate: coset DOS reconciliation (complete)

The coset decomposition makes the 6×6-coset downfolded problem strictly equivalent to a periodic defect array = the
SOC supercell. Reconciled directly against socsup's supercell eigenvalues (1040 each; O$_\mathrm{S}$ VBM = state 936,
V$_\mathrm{S}$ state 930):

**Three-way overview** (black = ground truth, blue = 10 bands, red = 22 bands; the V$_\mathrm{S}$ inset is a gap
zoom — the 10-band spurious peaks and misplaced doublet versus the 22-band clean fit settle it in one figure):

![3-way DOS comparison](../assets/socdos_3way.png)

![OS DOS truth 22b](../assets/socdos_OS_22b.png)

![VS DOS truth 22b](../assets/socdos_VS_22b.png)

**Manifold-width arbitration** (same host, same edmat; the only variable is the active window):

| V$_\mathrm{S}$ deep in-gap doublet | 10 bands | **22 bands** | supercell truth |
|---|---|---|---|
| position (eV) | ~0.93/1.0 (off by ~150 meV) | **+1.118/+1.166** | 1.086/~1.15 |
| spurious in-gap peaks | a string of false peaks, 0.18–0.73 | **none** | none |
| O$_\mathrm{S}$ CB onset | +1.538 (28 meV off) | +1.554 (12 meV off) | +1.566 |

**Conclusion and division of labour**: the V$_\mathrm{S}$ vacancy's deep state has strong S-p character, and the
Mo-d-dominated 10-band manifold cannot span it (the non-SOC lesson of "5 bands loses to 11 bands" repeats under SOC).
The 10-band frame is fine for **band-edge K-valley zooms** (the §1 figures on this page are qualitatively reliable
and 30× faster); **deep-gap or ground-truth-level quantitative work requires 22 bands**. The 22-band residual of
~20–40 meV comes from the N$_S$=16 chain truncation and can be deepened on demand.

### 4b. 22-band 144k dilute limit (three defects, new)

A 12×12 (144 k) × 22-spinor-band chain (SVD compression `svd_tol=1e-4`, ranks 1032/1032/1080,
4.6–4.7 h/chain on 4 Banff nodes) gives the dilute-limit DOS at $n_d=1/144$, compared three ways against the
supercell ground truth and the 6×6 array (η=10 meV):

![VS 22b 3-way](../assets/socdos22_VS.png)

![OS 22b 3-way](../assets/socdos22_OS.png)

![SES 22b 3-way](../assets/socdos22_SES.png)

| observable | truth (array 1/36) | coset 6×6 22b | **dilute 144k 22b** |
|---|---|---|---|
| V$_\mathrm{S}$ deep in-gap doublet (eV) | +1.085 / +1.135 | +1.115 / +1.165 | **+1.120 / +1.170** |
| O$_\mathrm{S}$ CB resonance (eV) | +1.545 | +1.545 | **+1.545** |
| Se$_\mathrm{S}$ gap | clean | — | **clean** (first downfolded DOS) |

**Conclusions**: (i) the dilution shift of the V$_\mathrm{S}$ deep state is ≤5 meV — **the deep in-gap levels have
already reached the isolated-defect limit at 1/36 concentration** (defect–defect coupling is negligible); the
~30 meV array-vs-truth residual comes from the $N_S=16$ chain truncation and is concentration-independent.
(ii) O$_\mathrm{S}$'s CB resonance coincides across all three, cleanly concentration-independent.
(iii) The dilute spectrum uniquely resolves V$_\mathrm{S}$'s shallow +0.19 eV in-gap feature (the 22-band
confirmation of the +0.181 seen in the 10-band era). The chain host-VBM (−5.8859 eV) and the spectral-function
pipeline's GATE_VBM verify each other to the digit.

## 5. To do

- ~~Se$_\mathrm{S}$ SOC chain + ground-truth reconciliation + three-defect zooms~~ (all complete: the 22-band three-defect version in §1, DOS in §4b);
- ~~22-band cross-check of the band-edge zooms~~ (complete: the 22-band version in §1 is exactly that; V$_\mathrm{S}$'s CB flat bands go from 10-band +1.694/+1.758 to a single 22-band feature at +1.698, Se$_\mathrm{S}$'s CB is clean, and on the V$_\mathrm{S}$ VBM side −0.213/−0.181/−0.036 plus a +0.185 in-gap spike (mutually confirmed by the +0.19 of the §4b dilute DOS));
- leave-one-out interpolation error for the 10-band frame (the Wannier-quality gate).
