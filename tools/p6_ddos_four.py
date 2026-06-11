#!/usr/bin/env python3
"""Defect-induced ΔDOS, FOUR treatments on the common grid of p6_ddos_three:
  (1) explicit-21 bare M   (2) 11-band bare M   (3) 11-band + 2nd-order Σ (over-screened)
  (4) 11-band + FULL-ORDER Σ (direct resolvent from explicit-60: vtilde_block_fesh60.dat)
Panels 1-3 reloaded from p6_ddos_three.npz; only (4) is computed here."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05

d=np.load("p6_ddos_three.npz")
oms=d["oms"]; dd={1:d["ddos1"],2:d["ddos2"],3:d["ddos3"]}
VBM=float(d["VBM"]); CBM=float(d["CBM"])

fb=FortranFile("vtilde_block_fesh60.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); fb.read_reals(np.float64); fb.read_ints(np.int32)
eps=fb.read_reals(np.float64)*RY
fb.read_record(np.complex128); fb.read_record(np.complex128)
Vt=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
print(f"fesh60 block: N_A={N_A}",flush=True)

V=(1.0/nk)*Vt*RY; I=np.eye(N_A); out=np.empty(len(oms))
for j,om in enumerate(oms):
    GA=1.0/(om+1j*ETA-eps)
    T=np.linalg.solve(I-V*GA[None,:],V)
    out[j]=-(1/np.pi)*np.imag(np.sum(GA*GA*np.diag(T)))
    if j%60==0: print(f"  omega {j}/{len(oms)}",flush=True)
dd[4]=out
np.savez("p6_ddos_four.npz",oms=oms,ddos1=dd[1],ddos2=dd[2],ddos3=dd[3],ddos4=dd[4],VBM=VBM,CBM=CBM)

def peak(y):
    m=(oms>VBM+0.05)&(oms<CBM-0.02); i=np.argmax(np.where(m,y,-1e9)); return oms[i],y[i]
cases=[(1,"(1) explicit summation — 21 bands, bare $M$","#1f3b73"),
       (2,"(2) active space — 11 bands, bare $M$","#2a7a2a"),
       (3,"(3) 11 bands + 2nd-order rest $\\Sigma^{(2)}$ — over-screened","#c1121f"),
       (4,"(4) 11 bands + FULL-ORDER rest $\\Sigma_{\\rm full}$ (direct resolvent)","#6a4c93")]
fig,axs=plt.subplots(4,1,figsize=(8.6,11.4),sharex=True)
for ax,(i,ttl,c) in zip(axs,cases):
    y=dd[i]
    ax.axvspan(oms[0],VBM,color="#cfe0f5",alpha=.5); ax.axvspan(CBM,oms[-1],color="#f6d3cf",alpha=.5)
    ax.axhline(0,color="#bbb",lw=.6); ax.plot(oms,y,c=c,lw=1.8)
    ax.fill_between(oms,0,y,where=(y>0),color=c,alpha=.18)
    pe,ph=peak(y)
    ax.annotate(f"$e$ peak @ {pe:.2f} eV",xy=(pe,ph),xytext=(pe+0.13,max(ph*0.6,6)),
                fontsize=9,color=c,arrowprops=dict(arrowstyle="->",color=c))
    for e,nm,col in [(-5.83,"DFT $a_1$","#666"),(-4.77,"DFT $e$","#999")]:
        ax.axvline(e,color=col,lw=.9,ls=":")
        if i==1: ax.text(e,ax.get_ylim()[1],nm,color=col,fontsize=8,rotation=90,va="top",ha="right")
    ax.set_title(ttl,fontsize=10.5,loc="left"); ax.set_ylabel("$\\Delta$DOS (1/eV)")
axs[-1].set_xlim(VBM-0.25,CBM+0.15); axs[-1].set_xlabel("$\\omega$ (eV)")
fig.suptitle("Defect-induced $\\Delta$DOS — the rest-dressing story in four panels ($\\eta=50$ meV)",fontsize=12,y=0.995)
plt.tight_layout(); plt.savefig("p6_ddos_four.png",dpi=140,bbox_inches="tight")
print("wrote p6_ddos_four.png")
for i,t,_ in cases:
    pe,ph=peak(dd[i]); print(f"panel {i}: e peak {pe:.3f} eV (VBM+{pe-VBM:.2f}), height {ph:.1f}/eV")
