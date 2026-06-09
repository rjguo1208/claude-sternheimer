#!/usr/bin/env python3
"""Publication figure for the explicit (non-Wannier) 21-band T-matrix: (left) DOS in the gap with the
in-gap e peak, (right) defect-level diagram explicit-T-matrix vs DFT supercell.  Replots from the npz."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("p6_explicit_dos.npz")
oms=d["oms"]; dos=d["dos"]; dos0=d["dos0"]; VBM=float(d["VBM"]); CBM=float(d["CBM"])
ea1,ee=-5.930,-4.588                       # explicit T-matrix levels (eig of E+M/N_k)
da1,de=-5.83,-4.77                         # DFT supercell levels
fig,(axL,axR)=plt.subplots(1,2,figsize=(11,4.8),gridspec_kw={'width_ratios':[2.1,1]})

# ---- left: DOS (zoomed on the gap so the e peak is visible) ----
axL.axvspan(oms[0],VBM,color="#cfe0f5",alpha=.55); axL.axvspan(CBM,oms[-1],color="#f6d3cf",alpha=.55)
axL.plot(oms,dos,c="#1f3b73",lw=1.9,label="DOS with defect (explicit 21-band T-matrix)")
axL.plot(oms,dos0,c="#888",lw=1.1,ls="--",label="bare host DOS (no defect)")
axL.axvline(ee,c="#e36414",lw=1.1,ls=":"); axL.axvline(de,c="#e36414",lw=1.0,ls=(0,(1,2)),alpha=.55)
axL.annotate("$e$ defect level\n(in-gap peak)",xy=(ee,12.6),xytext=(ee-0.05,42),ha="center",fontsize=9,
             color="#e36414",arrowprops=dict(arrowstyle="->",color="#e36414"))
axL.text(VBM-0.05,52,"valence",ha="right",fontsize=8,color="#26408b")
axL.text(CBM+0.05,52,"conduction",ha="left",fontsize=8,color="#9b2226")
axL.set_ylim(0,58); axL.set_xlim(VBM-0.8,CBM+0.8)
axL.set_xlabel("$\\omega$ (eV)"); axL.set_ylabel("DOS (states/eV)")
axL.set_title("Explicit non-Wannier T-matrix DOS (gap, zoomed)",fontsize=10)
axL.legend(fontsize=8,loc="upper center")

# ---- right: level diagram, explicit vs DFT ----
lo,hi=-6.25,-4.05; xA=(0.08,0.42); xB=(0.58,0.92)
for (x0,x1),ttl in [(xA,"explicit\nT-matrix"),(xB,"DFT\nsupercell")]:
    axR.axhspan(lo,VBM,xmin=x0,xmax=x1,color="#cfe0f5",alpha=.55)
    axR.axhspan(CBM,hi,xmin=x0,xmax=x1,color="#f6d3cf",alpha=.55)
    axR.text((x0+x1)/2,hi-0.05,ttl,ha="center",va="top",fontsize=9,weight="bold")
for x,a1,e,tag in [(xA,ea1,ee,"expl"),(xB,da1,de,"dft")]:
    axR.hlines(a1,x[0],x[1],color="#c1121f",lw=2.6)
    axR.hlines(e,x[0],x[1],color="#e36414",lw=2.6)
    axR.hlines(e+0.018,x[0],x[1],color="#e36414",lw=2.6)          # doublet (two lines)
axR.text(xA[1]+0.01,ea1,f"$a_1$ {ea1:.2f}",va="center",fontsize=8,color="#c1121f")
axR.text(xA[1]+0.01,ee,f"$e$×2 {ee:.2f}",va="center",fontsize=8,color="#e36414")
axR.text(xB[1]+0.01,da1,f"{da1:.2f}",va="center",fontsize=8,color="#c1121f")
axR.text(xB[1]+0.01,de,f"{de:.2f}",va="center",fontsize=8,color="#e36414")
axR.plot([xA[1],xB[0]],[ea1,da1],c="#bbb",lw=.7,ls=(0,(1,2))); axR.plot([xA[1],xB[0]],[ee,de],c="#bbb",lw=.7,ls=(0,(1,2)))
axR.set_xlim(0,1.15); axR.set_ylim(lo,hi); axR.set_xticks([]); axR.set_ylabel("Energy (eV)")
axR.set_title("Defect levels: $a_1$ + $e$ doublet",fontsize=10)
plt.tight_layout(); plt.savefig("p6_explicit_dos_fig.png",dpi=140,bbox_inches="tight"); print("wrote p6_explicit_dos_fig.png")
