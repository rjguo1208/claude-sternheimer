#!/usr/bin/env python3
"""DOS by DIRECT diagonalization of the single-defect Hamiltonian on the 12x12 grid:
  H_eff = diag(eps_i) + g,   g_ij = <i|dV|j> = M_ij / N_k   (M = stored raw block)
eigvalsh(H_eff) -> eigenvalues E_alpha -> DOS(w) = sum_alpha Lorentzian(w-E_alpha, eta).
Bare host DOS from eps_i for reference. Done for vacancy / O_S / Se_S (explicit-60,
N_A=8758). The in-gap eigenvalues are the bound states (delta-like peaks in the gap)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05; MID=-5.1
RUN="/anvil/scratch/x-rg47749/claude/sternheimerED/edt/run"
CASES=[("vtilde_explicit60.dat","S vacancy","#1f3b73"),
       ("vtilde_explicit60_OS.dat","O$_S$ substitution","#c1121f"),
       ("vtilde_explicit60_SeS.dat","Se$_S$ substitution","#2a7a2a")]

def dos_of(fn):
    fb=FortranFile(f"{RUN}/{fn}","r")
    NA,nk,_,_=fb.read_ints(np.int32); fb.read_reals(np.float64)
    idx=fb.read_ints(np.int32); a2b=idx[NA:]; eact=fb.read_reals(np.float64)
    M=fb.read_record(np.complex128).reshape((NA,NA),order="F"); fb.close()
    eps=eact*RY; Nk=float(nk)
    H=np.diag(eps.astype(complex))+M*RY/Nk; H=(H+H.conj().T)/2
    print(f"  {fn}: diagonalizing {NA}^2 ...",flush=True)
    w=np.linalg.eigvalsh(H)
    VBM=eps[(a2b<=13)&(eps<MID)].max(); CBM=eps[(a2b>=14)&(eps>MID)].min()
    ig=w[(w>VBM+3e-3)&(w<CBM-3e-3)]
    return w,eps,VBM,CBM,ig,Nk

oms=np.linspace(-7.5,-2.5,1000)
def lor(centers,Nk): return (ETA/np.pi/((oms[:,None]-centers[None,:])**2+ETA**2)).sum(1)/Nk

fig,axs=plt.subplots(3,1,figsize=(8.5,10),sharex=True)
res={}
for ax,(fn,lab,c) in zip(axs,CASES):
    w,eps,VBM,CBM,ig,Nk=dos_of(fn)
    dos=lor(w,Nk); dos0=lor(eps,Nk); res[lab]=(VBM,CBM,ig)
    ax.axvspan(VBM,CBM,color="#eee",alpha=.8,zorder=0)
    ax.plot(oms,dos0,c="#999",lw=1.0,ls="--",label="bare host (diag $H_0$)")
    ax.plot(oms,dos,c=c,lw=1.8,label=f"{lab} (diag $H_0+g$)")
    for e in ig: ax.axvline(e,c=c,lw=0.8,ls=":",alpha=.7)
    for E,nm in [(-5.83,"DFT $a_1$"),(-4.77,"DFT $e$")]:
        ax.axvline(E,c="#444",lw=0.8,ls=":",alpha=.5)
    ax.set_title(f"{lab}   (VBM {VBM:.2f}, CBM {CBM:.2f}; in-gap: "
                 +(", ".join(f"{x-VBM:+.2f}" for x in ig) if len(ig) else "none")+")",
                 fontsize=10,loc="left")
    ax.set_ylabel("DOS (states/eV/cell)"); ax.legend(fontsize=8,loc="upper right")
    ax.set_ylim(0,None)
axs[-1].set_xlim(-7.5,-2.5); axs[-1].set_xlabel("$\\omega$ (eV)")
fig.suptitle("DOS by direct diagonalization of $H_{\\rm eff}=\\mathrm{diag}(\\varepsilon)+g$  "
             "(12$\\times$12, explicit-60; gap shaded; $\\eta$=50 meV)",y=0.995)
plt.tight_layout(); plt.savefig(f"{RUN}/p7_dos_diag.png",dpi=140,bbox_inches="tight")
print("\nwrote p7_dos_diag.png")
for lab,(VBM,CBM,ig) in res.items():
    print(f"  {lab}: gap [{VBM:.3f},{CBM:.3f}]  in-gap levels (rel VBM): "
          +(", ".join(f"{x-VBM:+.3f}" for x in ig) if len(ig) else "NONE"))
