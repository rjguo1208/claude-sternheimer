#!/usr/bin/env python3
"""Defect-induced ΔDOS(ω) = -1/π Im Tr[G^A T G^A] for three treatments, on one common (ω-grid, η):
  (1) explicit summation, 21 bands, raw Bloch (bare M, no rest dressing)   [from p6_explicit_dos.npz]
  (2) active-space 11 bands (bare M, no rest dressing)
  (3) active-space 11 bands + 2nd-order rest space (dressed Ṽ = M + Σ, Σ from the full bands)
ΔDOS ≈ 0 in the gap for the bare host, so the in-gap features ARE the defect levels (a1, e)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05

# --- case 1: explicit 21-band (already computed, same grid & eta) ---
d=np.load("p6_explicit_dos.npz")
oms=d["oms"]; VBM=float(d["VBM"]); CBM=float(d["CBM"])
ddos1=d["dos"]-d["dos0"]

# --- load the 11-band block (ε, bare M, dressed Ṽ) ---
def load(fn):
    fb=FortranFile(fn,"r")
    N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32)
    fb.read_reals(np.float64); fb.read_ints(np.int32)
    eps=fb.read_reals(np.float64)*RY; N=int(N_A)
    M=fb.read_record(np.complex128).reshape((N,N),order="F")
    fb.read_record(np.complex128)                       # Σ
    Vt=fb.read_record(np.complex128).reshape((N,N),order="F")
    fb.close(); return int(nk),eps,M,Vt
nk,eps,M,Vt=load("vtilde_block.dat")

def ddos(eps,nk,mat,oms,eta):
    N=len(eps); I=np.eye(N); V=(1.0/nk)*mat*RY; out=np.empty(len(oms))
    for j,om in enumerate(oms):
        GA=1.0/(om+1j*eta-eps)
        T=np.linalg.solve(I-V*GA[None,:],V)             # T=[1-VG^A]^{-1}V
        out[j]=-(1/np.pi)*np.imag(np.sum(GA*GA*np.diag(T)))   # Tr[G^A T G^A]
    return out
ddos2=ddos(eps,nk,M,oms,ETA)                            # 11-band bare M
ddos3=ddos(eps,nk,Vt,oms,ETA)                           # 11-band + 2nd-order rest
np.savez("p6_ddos_three.npz",oms=oms,ddos1=ddos1,ddos2=ddos2,ddos3=ddos3,VBM=VBM,CBM=CBM)

def peak(dd,lo=0.05,hi=0.02):
    m=(oms>VBM+lo)&(oms<CBM-hi); i=np.argmax(np.where(m,dd,-1e9)); return oms[i],dd[i]

cases=[("(1) explicit summation — 21 bands, bare $M$",ddos1,"#1f3b73"),
       ("(2) active space — 11 bands, bare $M$",ddos2,"#2a7a2a"),
       ("(3) 11 bands + 2nd-order rest space — $\\tilde V=M+\\Sigma$",ddos3,"#c1121f")]
fig,axs=plt.subplots(3,1,figsize=(8.6,9),sharex=True)
for ax,(ttl,dd,c) in zip(axs,cases):
    ax.axvspan(oms[0],VBM,color="#cfe0f5",alpha=.5); ax.axvspan(CBM,oms[-1],color="#f6d3cf",alpha=.5)
    ax.axhline(0,color="#bbb",lw=.6); ax.plot(oms,dd,c=c,lw=1.8)
    ax.fill_between(oms,0,dd,where=(dd>0),color=c,alpha=.18)
    pe,ph=peak(dd)
    ax.annotate(f"$e$ peak @ {pe:.2f} eV",xy=(pe,ph),xytext=(pe+0.12,ph*0.65),
                fontsize=9,color=c,arrowprops=dict(arrowstyle="->",color=c))
    for e,nm,col in [(-5.83,"DFT $a_1$","#666"),(-4.77,"DFT $e$","#999")]:
        ax.axvline(e,color=col,lw=.9,ls=":")
        if ax is axs[0]: ax.text(e,ax.get_ylim()[1],nm,color=col,fontsize=8,rotation=90,va="top",ha="right")
    ax.set_title(ttl,fontsize=10.5,loc="left"); ax.set_ylabel("$\\Delta$DOS (1/eV)")
    ax.text(VBM-0.02,ax.get_ylim()[1]*0.9,"VBM",ha="right",fontsize=7,color="#26408b")
    ax.text(CBM+0.02,ax.get_ylim()[1]*0.9,"CBM",ha="left",fontsize=7,color="#9b2226")
axs[-1].set_xlim(VBM-0.25,CBM+0.15); axs[-1].set_xlabel("$\\omega$ (eV)")
fig.suptitle("Defect-induced $\\Delta$DOS of the MoS$_2$ S-vacancy (12$\\times$12, $\\eta=50$ meV)",fontsize=12,y=0.997)
plt.tight_layout(); plt.savefig("p6_ddos_three.png",dpi=140,bbox_inches="tight")
print("wrote p6_ddos_three.png")
print(f"e-peak  (1) explicit : {peak(ddos1)[0]:.3f} eV  (VBM+{peak(ddos1)[0]-VBM:.2f})")
print(f"e-peak  (2) bare-M 11: {peak(ddos2)[0]:.3f} eV  (VBM+{peak(ddos2)[0]-VBM:.2f})")
print(f"e-peak  (3) dressed  : {peak(ddos3)[0]:.3f} eV  (VBM+{peak(ddos3)[0]-VBM:.2f})")
