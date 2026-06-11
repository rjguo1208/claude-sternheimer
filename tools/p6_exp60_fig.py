#!/usr/bin/env python3
"""Publication figure: (left) band convergence of the explicit e level onto DFT;
(right) explicit-60 DOS zoomed on the gap with the in-gap levels.  Reads
p6_explicit60_levels.npz (run dir)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("p6_explicit60_levels.npz")
oms=d["oms"]; dos=d["dos"]; dos0=d["dos0"]; VBM=float(d["VBM"]); CBM=float(d["CBM"]); ig=d["ig"]
nb=[11,21,60.8]; ee=[1.4949,1.3482,1.2051]; dft=1.1856
fig,(axL,axR)=plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={'width_ratios':[1,1.4]})

# left: convergence curve
axL.plot(nb,ee,"o-",c="#1f3b73",lw=1.8,ms=7)
axL.axhline(dft,c="#e36414",lw=1.4,ls="--")
axL.text(12,dft-0.02,"DFT supercell  +1.19",color="#e36414",fontsize=9,va="top")
for x,y,lab in zip(nb,ee,["11","21","60"]):
    axL.annotate(f"+{y:.3f}",(x,y),textcoords="offset points",xytext=(8,6),fontsize=9,color="#1f3b73")
axL.set_xlabel("explicit bands per $k$"); axL.set_ylabel("$e$ level above VBM (eV)")
axL.set_title("Band convergence of the explicit $e$ level",fontsize=10)
axL.set_xlim(5,70); axL.set_ylim(1.12,1.56); axL.grid(alpha=.3)

# right: explicit-60 gap DOS
axR.axvspan(oms[0],VBM,color="#cfe0f5",alpha=.55); axR.axvspan(CBM,oms[-1],color="#f6d3cf",alpha=.55)
axR.plot(oms,dos,c="#1f3b73",lw=1.8,label="explicit-60 ($N_A$=8758)")
axR.plot(oms,dos0,c="#888",lw=1.0,ls="--",label="bare host")
for e,nm,c in [(-5.83,"DFT $a_1$","#c1121f"),(-4.77,"DFT $e$","#e36414")]:
    axR.axvline(e,c=c,lw=1.2,ls=":")
    axR.text(e,55,nm,color=c,fontsize=8,rotation=90,va="top",ha="right")
for e in ig: axR.axvline(e,c="#2a7a2a",lw=1.2,ls="--",alpha=.8)
axR.annotate("$e$ doublet $-4.731$\n(DFT $-4.77$)",xy=(float(ig.max()),18),xytext=(-5.35,40),
             fontsize=9,color="#2a7a2a",arrowprops=dict(arrowstyle="->",color="#2a7a2a"))
axR.set_xlim(VBM-0.7,CBM+0.7); axR.set_ylim(0,60)
axR.set_xlabel("$\\omega$ (eV)"); axR.set_ylabel("DOS (states/eV)")
axR.set_title("Explicit-60 DOS (gap, zoomed): levels (green) vs DFT (dotted)",fontsize=10)
axR.legend(fontsize=8,loc="upper center")
plt.tight_layout(); plt.savefig("p6_exp60_fig.png",dpi=140,bbox_inches="tight")
print("wrote p6_exp60_fig.png")
