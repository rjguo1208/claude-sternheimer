#!/usr/bin/env python3
"""Explicit-60 defect levels + DOS via ONE eigendecomposition (no T-matrix solves).
H_eff = diag(eps) + M/N_k; in-gap eigenvalues = defect levels; DOS = Lorentzian sum
over eigenvalues (equivalent to the T-matrix DOS for the same H_eff).  Prints the
band-convergence line: e(11)=+1.50, e(21)=+1.35, e(60)=?  (DFT +1.19)."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05; MID=-5.1
FN=sys.argv[1] if len(sys.argv)>1 else "vtilde_explicit60.dat"

fb=FortranFile(FN,"r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); om0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:2*N_A]
eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
eps=eact*RY; Nk=float(nk)
print(f"N_A={N_A} nk={nk} bands/k~{N_A/nk:.1f}  NSCF band range {a2band.min()}-{a2band.max()}")
VBM=eps[(a2band<=13)&(eps<MID)].max(); CBM=eps[(a2band>=14)&(eps>MID)].min()
print(f"VBM={VBM:.3f}  CBM={CBM:.3f}  gap={CBM-VBM:.3f} eV   (DFT ref: a1 -5.83/+0.13, e -4.77/+1.19)")

H=np.diag(eps.astype(complex))+(1.0/Nk)*M*RY; H=(H+H.conj().T)/2
print("diagonalizing H_eff ...", flush=True)
w=np.linalg.eigvalsh(H)
ig=w[(w>VBM+3e-3)&(w<CBM-3e-3)]
print(f"in-gap levels ({len(ig)}):")
for e in ig:
    # degeneracy within 1e-3 eV
    deg=int(np.sum(np.abs(w-e)<1e-3))
    print(f"   E = {e:.4f} eV   (+{e-VBM:.3f} above VBM)   multiplicity~{deg}")
print(f"\nband convergence of e:  11bnd +1.50  ->  21bnd +1.35  ->  60bnd +{ig.max()-VBM:.3f}   (DFT +1.19)")

# DOS from eigenvalues (Lorentzian sum) vs bare host
WLO,WHI=VBM-1.0,CBM+1.0; oms=np.linspace(WLO,WHI,500)
dos =np.sum(ETA/np.pi/((oms[:,None]-w[None,:])**2+ETA**2),axis=1)
dos0=np.sum(ETA/np.pi/((oms[:,None]-eps[None,:])**2+ETA**2),axis=1)
np.savez("p6_explicit60_levels.npz",w=w,eps=eps,oms=oms,dos=dos,dos0=dos0,VBM=VBM,CBM=CBM,ig=ig)

fig,ax=plt.subplots(figsize=(8,5))
ax.axvspan(WLO,VBM,color="#cfe0f5",alpha=.5); ax.axvspan(CBM,WHI,color="#f6d3cf",alpha=.5)
ax.plot(oms,dos,c="#1f3b73",lw=1.6,label="explicit-60 H_eff DOS")
ax.plot(oms,dos0,c="#999",lw=1.0,ls="--",label="bare host")
for e,nm,c in [(-5.83,"DFT $a_1$","#c1121f"),(-4.77,"DFT $e$","#e36414")]:
    ax.axvline(e,c=c,lw=1.2,ls=":")
for e in ig: ax.axvline(e,c="#2a7a2a",lw=1.0,ls="--",alpha=.7)
ax.set_xlim(WLO,WHI); ax.set_xlabel("$\\omega$ (eV)"); ax.set_ylabel("DOS (states/eV)")
ax.set_title(f"Explicit-60 (N_A={N_A}): in-gap levels (green) vs DFT (dotted)")
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig("p6_explicit60_levels.png",dpi=140,bbox_inches="tight")
print("wrote p6_explicit60_levels.png")
