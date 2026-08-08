#!/usr/bin/env python3
"""Does the defect vertex survive Wannier interpolation?

(a) Locality of the pair kernel M(R_e,R_p) about the defect -- the assumption the
    whole interpolation rests on.  Built from the 6x6 sub-grid and from the full
    12x12 set of the same 144-k run, so the two curves are the same quantity
    measured with two BvK cells.
(b) Leave-one-out test: build the kernel from the Gamma coset alone (6.2% of the
    pairs), predict all 144x144, compare element-wise in the same Wannier gauge."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RY = 13.605693122994; NK = 144
z = np.load("/home/rjguo/edi_tmatrix/sternA/vsrlx/loo.npz")
C12, C66, CG = "#1f3b73", "#c1121f", "#8a8f98"
fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.3))

# ---- (a) locality ------------------------------------------------------------
def shell(d, v, edges):
    """envelope: the largest element in each shell, shells centred on integers so
    the on-site term R_e = R_def is its own point and not lumped with |R| = 1."""
    x, y = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (d >= lo) & (d < hi)
        if m.any(): x.append(0.5*(lo+hi) if lo > 0 else 0.0); y.append(v[m].max())
    return np.array(x), np.array(y)
ed = np.arange(-0.5, 8.6, 1.0)
for k, lab, c, ls, mk in (("66", "6$\\times$6 BvK cell", C66, "--", "s"),
                          ("1212", "12$\\times$12 BvK cell", C12, "-", "o")):
    x, y = shell(z["d_"+k], z["max_"+k], ed)
    ax[0].semilogy(x, y, ls, color=c, lw=1.9, marker=mk, ms=5.5, label=lab)
for r, c in ((3.0, C66), (6.0, C12)):
    ax[0].axvline(r, color=c, lw=.9, ls=":", alpha=.7)
ax[0].text(3.0, 2.0, " 6$\\times$6 cell edge", color=C66, fontsize=8, ha="left")
ax[0].text(6.0, 2.0, " 12$\\times$12 edge", color=C12, fontsize=8, ha="left")
ax[0].annotate("$R_e=R_\\mathrm{def}$", xy=(0.0, 1.46), xytext=(0.9, 2.6), fontsize=8.5,
               color="#444", arrowprops=dict(arrowstyle="->", color="#444", lw=.8))
ax[0].set_xlabel("$|R_e - R_\\mathrm{def}|$   (lattice constants)")
ax[0].set_ylabel("$\\max_{R_p}\\,|\\mathcal{M}(R_e,R_p)|$")
ax[0].set_xlim(-0.4, 7.6); ax[0].set_ylim(1e-4, 5.0)
ax[0].legend(frameon=False, fontsize=9)
ax[0].set_title("(a)  locality of the pair kernel about the defect", fontsize=10, loc="left")

# ---- (b) leave-one-out -------------------------------------------------------
err, on = z["err"], z["on"]
mev = 1e3*RY*err/NK
grp = (("one off-grid", np.outer(on, ~on) | np.outer(~on, on), C66, "--"),
       ("both off-grid (93.8% of pairs)", np.outer(~on, ~on), C12, "-"))
for lab, m, c, ls in grp:
    v = np.sort(mev[m].ravel()); v = np.maximum(v, 1e-9)
    ax[1].semilogx(v, np.linspace(0, 100, len(v)), ls, color=c, lw=2.0, label=lab)
for x, lab, c in ((7.6, "downfold residual", "#7a5c00"), (11.0, "6$\\times$6 image error", "#2a7a2a")):
    ax[1].axvline(x, color=c, lw=1.1, ls="-.", alpha=.85)
    ax[1].text(x*1.09, 28, lab, color=c, fontsize=8.5, rotation=90, va="center", ha="left")
ax[1].text(0.13, 92, "training pairs reproduce to $4\\times10^{-15}$\n(off scale to the left)",
           fontsize=8.5, color="#666")
ax[1].set_xlim(0.12, 40); ax[1].set_ylim(0, 100)
ax[1].set_xlabel("interpolation error of one vertex element   (meV in $H_\\mathrm{eff}$)")
ax[1].set_ylabel("percentile of $(k,k')$ pairs")
ax[1].legend(frameon=False, fontsize=9, loc="center left", bbox_to_anchor=(0.02, 0.42))
ax[1].set_title("(b)  leave-one-out: predict the 93.8% never seen", fontsize=10, loc="left")
for a in ax: a.tick_params(labelsize=9); a.grid(alpha=.18, lw=.6)
fig.tight_layout()
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/wannier_decay.png"
fig.savefig(out, dpi=170); print("wrote", out)
for lab, m, c, ls in grp:
    v = mev[m]
    print("  %-38s median %7.3f  p90 %7.3f  max %7.3f meV"
          % (lab.replace("$\\times$", "x"), np.median(v), np.percentile(v, 90), v.max()))
