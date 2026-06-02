#!/usr/bin/env python3
"""Intraband (top valence band) scattering map: |M(k_f,K)| and |Vtilde(k_f,K)| from the
high-symmetry K point to all 144 final k_f, in Cartesian k (folded to the 1st BZ)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:]; fb.read_reals(np.float64)
M =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
Sg=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
V =fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F"); fw.close()   # crystal coords

# reciprocal lattice (from the primitive lattice = V_d.cube supercell / 6), 1/Bohr
a1=np.array([240*0.150478, 0.0])/6;  a2=np.array([240*-0.075239, 240*0.130318])/6
area=a1[0]*a2[1]-a1[1]*a2[0]
b1=2*np.pi/area*np.array([ a2[1],-a2[0]]); b2=2*np.pi/area*np.array([-a1[1], a1[0]])
def kcart_BZ(k1,k2):                     # fold to 1st BZ (nearest to Gamma)
    best=None
    for p in (-1,0,1):
        for q in (-1,0,1):
            c=(k1+p)*b1+(k2+q)*b2
            if best is None or c@c < best@best: best=c
    return best
kc=np.array([kcart_BZ(xk[0,k],xk[1,k]) for k in range(nk)])         # [k,2] Cartesian

# band 17 (top valence) a-index for each k ; K = (1/3,1/3)
b_top=17
a_of_k=np.full(nk,-1,int)
for a in range(N_A):
    if a2band[a]==b_top: a_of_k[a2k[a]-1]=a
Kidx=int(np.argmin([min(((xk[0,k]-1/3+p)**2+(xk[1,k]-1/3+q)**2) for p in(-1,0,1) for q in(-1,0,1)) for k in range(nk)]))
a_i=a_of_k[Kidx]
print(f"K at k-index {Kidx+1}, crystal=({xk[0,Kidx]:.3f},{xk[1,Kidx]:.3f}); band={b_top}")
Mv=np.abs(M[a_of_k,a_i]); Vv=np.abs(V[a_of_k,a_i])
print(f"|M(K,K)|={Mv[Kidx]:.3f}  |Vtilde(K,K)|={Vv[Kidx]:.3f} Ry ;  max|M|={Mv.max():.3f} max|Vtilde|={Vv.max():.3f}")

# high-symmetry points for overlay
G=np.array([0,0]); Kpt=kcart_BZ(1/3,1/3); Mpt=kcart_BZ(1/2,0)
fig,ax=plt.subplots(1,2,figsize=(12.5,5.2))
for p,(val,ttl) in enumerate([(Mv,r"$|M(k_f,K)|$  (Born)"),(Vv,r"$|\tilde V(k_f,K)|$  (downfolded)")]):
    sc=ax[p].scatter(kc[:,0],kc[:,1],c=val,s=210,marker="h",cmap="viridis",edgecolors="0.6",linewidths=0.3)
    plt.colorbar(sc,ax=ax[p],label="Ry")
    for pt,lab in [(G,"Γ"),(Kpt,"K (source)"),(Mpt,"M")]:
        ax[p].plot(*pt,"r*" if "source" in lab else "w+",ms=14 if "source" in lab else 9,mew=2)
        ax[p].annotate(lab,pt,textcoords="offset points",xytext=(6,6),color="r" if "source" in lab else "w",fontsize=9,fontweight="bold")
    ax[p].set_xlabel("$k_x$ (1/Bohr)"); ax[p].set_ylabel("$k_y$ (1/Bohr)")
    ax[p].set_title(ttl+f"\nintraband, top valence band (b={b_top})"); ax[p].set_aspect("equal"); ax[p].grid(alpha=.2)
plt.tight_layout(); plt.savefig("p5b_kmap.png",dpi=120); print("wrote p5b_kmap.png")
