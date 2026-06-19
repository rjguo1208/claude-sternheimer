#!/usr/bin/env python3
"""Band-convergence of the in-gap defect level from the SAME explicit-60 M (no new edt run):
slice the host-band basis to a2band<=thr, build H_eff=diag(eps)+M/N_k on the subset, eigvalsh,
track the in-gap eigenvalue(s) vs #bands. A REAL bound state converges (vacancy e -> DFT +1.19);
a basis-truncation artifact keeps descending toward the VBM and would leave the gap in a complete
basis (suspected for the isovalent O_S, whose DFT shows NO in-gap state)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; MID=-5.1
RUN="/anvil/scratch/x-rg47749/claude/sternheimerED/edt/run"
CASES=[("vtilde_explicit60.dat","S vacancy","#1f3b73"),
       ("vtilde_explicit60_OS.dat","O$_S$","#c1121f")]
THR=[17,19,22,27,32,42,52,67]    # a2band cutoffs -> ~11,13,16,21,26,36,46,61 bands

def load(fn):
    fb=FortranFile(f"{RUN}/{fn}","r"); NA,nk,_,_=fb.read_ints(np.int32); fb.read_reals(np.float64)
    idx=fb.read_ints(np.int32); a2b=idx[NA:]; eact=fb.read_reals(np.float64)
    M=fb.read_record(np.complex128).reshape((NA,NA),order="F"); fb.close()
    return nk,a2b,eact*RY,M

fig,ax=plt.subplots(figsize=(8,5)); res={}
for fn,lab,c in CASES:
    nk,a2b,eps,M=load(fn)
    VBM=eps[(a2b<=13)&(eps<MID)].max(); CBM=eps[(a2b>=14)&(eps>MID)].min()
    nbs=[]; tops=[]
    print(f"\n{lab} (VBM {VBM:.3f}, CBM {CBM:.3f}):")
    for thr in THR:
        sel=np.where(a2b<=thr)[0]; nb=len(sel)//nk
        Hs=np.diag(eps[sel].astype(complex))+M[np.ix_(sel,sel)]*RY/nk; Hs=(Hs+Hs.conj().T)/2
        w=np.linalg.eigvalsh(Hs); ig=w[(w>VBM+3e-3)&(w<CBM-3e-3)]
        deep=(ig.max()-VBM) if len(ig) else np.nan
        nbs.append(nb); tops.append(ig.max() if len(ig) else np.nan)
        print(f"  {nb:3d} bands: in-gap {len(ig):2d}, highest = {deep:+.3f} above VBM" if len(ig)
              else f"  {nb:3d} bands: NO in-gap state")
    res[lab]=(np.array(nbs),np.array(tops),VBM,CBM)
    ax.plot(nbs,np.array(tops)-VBM,'o-',c=c,label=lab)
ax.axhline(1.19,color="#1f3b73",ls=":",lw=1,label="DFT vacancy $e$ (+1.19)")
ax.axhline(0.0,color="#999",ls="--",lw=1,label="VBM")
ax.set_xlabel("# host bands in basis"); ax.set_ylabel("highest in-gap level $-$ VBM (eV)")
ax.set_title("Band convergence of the in-gap level (12$\\times$12, sliced from explicit-60 $M$)")
ax.legend(fontsize=9); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(f"{RUN}/p7_bandconv.png",dpi=140,bbox_inches="tight")
print("\nwrote p7_bandconv.png")
