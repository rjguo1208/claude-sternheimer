#!/usr/bin/env python3
"""Gamma-point energy-level diagram of the S-vacancy MoS2 supercell SCF, next to the active-space
T-matrix result, so the DFT defect levels (a1 singlet near VBM, e doublet in the gap) are visible and
directly comparable.  Parses eigenvalues + occupations from scf.out."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
F="/anvil/projects/x-che190065/rjguo/qe-7.5/edi-dev/test_edinterp/edi_run/defect_super/scf.out"
L=open(F).read().splitlines()
i0=max(i for i,l in enumerate(L) if "k = 0.0000 0.0000 0.0000" in l)
ev=[]; i=i0+1
while "occupation numbers" not in L[i]:
    ev+=[float(x) for x in L[i].split()]; i+=1
i+=1; occ=[]
while "Fermi" not in L[i] and "highest" not in L[i]:
    t=L[i].split()
    if t:
        try: occ+=[float(x) for x in t]
        except: pass
    i+=1
ev=np.array(ev); occ=np.array(occ[:len(ev)]); EF=-5.3012
homo=int(np.where(occ>0.5)[0][-1]); lumo=homo+1
VBM=ev[homo-1] if ev[homo]-ev[homo-1]>0.08 else ev[homo]
g=lumo
while g+1<len(ev) and ev[g+1]-ev[lumo]<0.03: g+=1
CBM=ev[g+1]; a1=ev[homo]; e_lvl=ev[lumo]; e_deg=g-lumo+1
print(f"HOMO state {homo+1} = {a1:.4f} eV (occ {occ[homo]:.2f})   LUMO state {lumo+1} = {e_lvl:.4f} eV x{e_deg}")
print(f"VBM(bulk)={VBM:.4f}  CBM(bulk)={CBM:.4f}  gap={CBM-VBM:.3f} eV   EF={EF}")
print(f"a1: {a1-VBM:+.3f} eV above VBM     e: {e_lvl-VBM:+.3f} above VBM ({e_lvl-CBM:+.3f} from CBM)")
for j in range(homo-5,lumo+6):
    print(f"   state {j+1}: {ev[j]:+.4f} eV  occ={occ[j]:.2f}")

# ----- active-space T-matrix result (this project) -----
VBM_t,CBM_t=-5.937,-4.279; res_t=-5.55     # with-rest defect resonance

fig,ax=plt.subplots(figsize=(8.2,7.6)); lo,hi=-6.4,-3.7
xA=(0.05,0.42); xB=(0.58,0.95)             # two columns
for (x0,x1),vb,cb,ttl in [(xA,VBM,CBM,"Supercell DFT\n(107-atom, Γ)"),
                          (xB,VBM_t,CBM_t,"Active-space\nT-matrix")]:
    ax.axhspan(lo,vb,xmin=x0,xmax=x1,color="#cfe0f5",alpha=.55,zorder=0)
    ax.axhspan(cb,hi,xmin=x0,xmax=x1,color="#f6d3cf",alpha=.55,zorder=0)
    ax.text((x0+x1)/2,hi-0.10,ttl,ha="center",va="top",fontsize=10,weight="bold",
            bbox=dict(facecolor="white",alpha=0.75,edgecolor="none",pad=1))
# supercell: manifold lines + defect levels
for j in range(len(ev)):
    if lo<=ev[j]<=hi and not (homo<=j<=g):
        ax.hlines(ev[j],*xA,color=("#26408b" if occ[j]>0.5 else "#9aa0a6"),lw=0.9,alpha=.85)
ax.hlines(a1,*xA,color="#c1121f",lw=3,zorder=6)
ax.hlines(e_lvl,xA[0],xA[1],color="#e36414",lw=3,zorder=6)
ax.text(xA[1]+0.01,a1,f"$a_1$ (occ.) {a1:.2f}",va="center",fontsize=9,color="#c1121f")
ax.text(xA[1]+0.01,e_lvl,f"$e$ (empty ×{e_deg}) {e_lvl:.2f}",va="bottom",fontsize=9,color="#e36414")
# T-matrix: the single resonance + a marker for the missing e
ax.hlines(res_t,*xB,color="#6a4c93",lw=3,zorder=6)
ax.text(xB[0]-0.01,res_t,f"{res_t:.2f}  ",va="center",ha="right",fontsize=9,color="#6a4c93")
ax.text((xB[0]+xB[1])/2,res_t-0.13,"resonance (with rest)",ha="center",fontsize=8,color="#6a4c93")
ax.text((xB[0]+xB[1])/2,e_lvl,"✗  $e$ absent",ha="center",va="center",fontsize=10,color="#e36414",
        bbox=dict(facecolor="white",alpha=0.7,edgecolor="#e36414",boxstyle="round,pad=0.2"))
# edges + Fermi
for (x0,x1),vb,cb in [(xA,VBM,CBM),(xB,VBM_t,CBM_t)]:
    ax.hlines(vb,x0,x1,color="#26408b",lw=1,ls=":"); ax.hlines(cb,x0,x1,color="#9b2226",lw=1,ls=":")
ax.hlines(EF,xA[0],xA[1],color="#2a9d8f",lw=1.4,ls="--")
ax.text(xA[0]+0.005,EF-0.10,f"$E_F$ {EF:.2f}",fontsize=8,color="#2a9d8f")
# guide line linking a1 <-> resonance
ax.plot([xA[1],xB[0]],[a1,res_t],color="#888",lw=0.8,ls=(0,(1,2)))
ax.text(0.5,lo+0.07,"VBM/CBM aligned (DFT −5.95/−4.24 ; T-mat −5.94/−4.28)",ha="center",fontsize=8,color="#555")
ax.set_xlim(0,1.18); ax.set_ylim(lo,hi); ax.set_xticks([]); ax.set_ylabel("Energy (eV)")
ax.set_title("Γ-point defect levels: S-vacancy MoS$_2$ — supercell DFT vs active-space T-matrix",fontsize=10)
plt.tight_layout(); plt.savefig("/anvil/scratch/x-rg47749/claude/sternheimerED/edt/run/supercell_levels.png",dpi=140,bbox_inches="tight")
print("wrote supercell_levels.png")
