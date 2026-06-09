#!/usr/bin/env python3
"""Explicit (non-Wannier) T-matrix + DOS from the raw 21-band Bloch matrix elements M=<nk|dV|mk'>
on the 12x12 grid.  Defect levels = in-gap eigenvalues of H_eff = diag(eps) + M/N_k; the DOS is
built from the T-matrix T=[1-V G^A]^{-1}V, V=M/N_k, G^A=diag(1/(w+ieta-eps))."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05    # eV
fb=FortranFile("vtilde_explicit20.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); om0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:2*N_A]
eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.read_record(np.complex128)                       # Sigma (=0, born_only)
fb.read_record(np.complex128); fb.close()           # Vtilde (=M)
eps=eact*RY; Nk=float(nk)
print(f"N_A={N_A} nk={nk} bands/k~{N_A/nk:.1f}  NSCF band range {a2band.min()}-{a2band.max()}")
val=a2band<=13; con=a2band>=14
VBM=eps[val].max(); CBM=eps[con].min()
print(f"VBM={VBM:.3f}  CBM={CBM:.3f}  gap={CBM-VBM:.3f} eV    (DFT ref: a1=-5.83 +0.13, e=-4.77 +1.19)")

# ---- defect levels = in-gap eigenvalues of H_eff = diag(eps) + c*M ----
for c,tag in [(1.0/Nk,"M/N_k (physical)"),(1.0,"M (no 1/N_k)")]:
    H=np.diag(eps.astype(complex))+c*M*RY; H=(H+H.conj().T)/2
    w=np.linalg.eigvalsh(H)
    ig=w[(w>VBM+2e-3)&(w<CBM-2e-3)]
    s=", ".join(f"{e:.3f}({e-VBM:+.2f}>VBM)" for e in ig) if len(ig) else "none"
    print(f"  H_eff[{tag}]: {len(ig):2d} in-gap level(s): {s}")

# ---- T-matrix DOS (V=M/N_k), gap-focused omega grid ----
V=(1.0/Nk)*M*RY; I=np.eye(N_A)
WLO,WHI=VBM-1.0,CBM+1.0; NW=300; oms=np.linspace(WLO,WHI,NW)
dos=np.empty(NW); dos0=np.empty(NW)
for j,om in enumerate(oms):
    GA=1.0/(om+1j*ETA-eps)
    T=np.linalg.solve(I-V*GA[None,:], V)            # T=[1-V G^A]^{-1} V
    dos[j] =-(1/np.pi)*np.imag(np.sum(GA)+np.sum(GA*GA*np.diag(T)))   # Tr[G^A+G^A T G^A]
    dos0[j]=-(1/np.pi)*np.imag(np.sum(GA))
np.savez("p6_explicit_dos.npz",oms=oms,dos=dos,dos0=dos0,VBM=VBM,CBM=CBM,eps=eps)
print(f"DOS done; in-gap (defect) DOS peak = {np.max((dos-dos0)[(oms>VBM+0.05)&(oms<CBM-0.05)]):.2f} /eV")

# ---- plot ----
fig,ax=plt.subplots(figsize=(8,5))
ax.axvspan(WLO,VBM,color="#cfe0f5",alpha=.5); ax.axvspan(CBM,WHI,color="#f6d3cf",alpha=.5)
ax.plot(oms,dos,c="#1f3b73",lw=1.6,label="DOS with defect (explicit 21-band T-matrix)")
ax.plot(oms,dos0,c="#999",lw=1.0,ls="--",label="bare host DOS (no defect)")
for e,nm,col in [(-5.83,"DFT $a_1$","#c1121f"),(-4.77,"DFT $e$","#e36414")]:
    ax.axvline(e,c=col,lw=1.4,ls=":"); ax.text(e,ax.get_ylim()[1]*0.9,nm,color=col,fontsize=9,rotation=90,va="top")
ax.axvline(VBM,c="#26408b",lw=0.8); ax.axvline(CBM,c="#9b2226",lw=0.8)
ax.set_xlim(WLO,WHI); ax.set_xlabel("$\\omega$ (eV)"); ax.set_ylabel("DOS (states/eV)")
ax.set_title("Explicit non-Wannier T-matrix DOS — S-vacancy, 21 bands, 12$\\times$12 (gap region)")
ax.legend(fontsize=8,loc="upper left")
plt.tight_layout(); plt.savefig("p6_explicit_dos.png",dpi=140,bbox_inches="tight"); print("wrote p6_explicit_dos.png")
