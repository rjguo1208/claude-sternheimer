#!/usr/bin/env python3
"""Wannier-interpolated (48x48) intraband scattering maps, fixed initial k = K, top valence band:
  |M(k_f,K)|        bare electron-defect matrix element (Born)
  |V~(k_f,K)|       partial T-matrix (rest-dressed downfolded potential = M + Sigma)
  |T_PP(k_f,K)|     full T-matrix [1 - V~ G^A]^{-1} V~  (on-shell, omega = eps_top(K))

H_W(k) is built from the SAME filukk U and the active band energies, H_W = U^dag diag(E) U, so the
band-projection eigenvectors V(k)=U^dag(k) are gauge-consistent with the matrix-element rotation
O^W=U^dag O U. The round-trip is then exact on the coarse grid (verified at K against the raw block);
on the fine grid it is the standard Marzari-Vanderbilt Wannier interpolation. |.| color maps in
Cartesian (BZ-folded) k."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
ETA=0.05/RY; RCUT=4; NF=48

fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]
eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.read_record(np.complex128)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5))
M4=M.reshape(nk,nb,nk,nb); V4=V.reshape(nk,nb,nk,nb); Ek=eact.reshape(nk,nb)   # [k,slot] Ry
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))                                        # [k, band, wann]

# --- coarse Wannier R-space M^W, V^W, and H_W from the SAME U (gauge-consistent) ---
ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n) for m in ms for n in ms]); nR=len(Rg)
F=np.exp(2j*np.pi*(Rg@xk[:2]))                                    # [R,k]
def toR(O4):
    OW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),O4,Ub,optimize=True)   # U^dag O U
    return np.einsum("Af,fpgw,Bg->ApBw",F,OW,F.conj(),optimize=True)/nk**2
MWR=toR(M4); VWR=toR(V4)
HWk=np.einsum("ksp,ks,ksw->kpw",Ub.conj(),Ek,Ub,optimize=True)               # H_W(k)=U^dag E U
HWR=np.einsum("Rk,kpw->Rpw",np.exp(-2j*np.pi*(Rg@xk[:2])),HWk,optimize=True)/nk
def Hk_grid(kf2): return np.einsum("kR,Rpw->kpw",np.exp(2j*np.pi*(kf2@Rg.T)),HWR,optimize=True)

Kc=np.array([1/3,1/3]); eK,vK=np.linalg.eigh(Hk_grid(Kc[None])[0])
top=6; omegaT=eK[top]; vKtop=vK[:,top]            # slot 6 = VBM = NSCF band 13
print(f"K=(1/3,1/3); eps_top(K)=omega={omegaT*RY:.4f} eV (on-shell); eta={ETA*RY:.3g} eV")

# --- T^W(R',R) via Koster-Slater at omega=omegaT (fine G^A), same H_W ---
Rd=int(np.argmax([np.linalg.norm(VWR[r,:,r,:]) for r in range(nR)]))
def msh(R): return ((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2, ((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
sel=[r for r in range(nR) if max(abs(msh(r)[0]),abs(msh(r)[1]))<=RCUT]; ns=len(sel); dim=nb*ns
Vsub=np.zeros((dim,dim),complex)
for a,Ra in enumerate(sel):
    for b,Rb in enumerate(sel): Vsub[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=VWR[Ra,:,Rb,:]
kfg=np.array([[i/NF,j/NF] for i in range(NF) for j in range(NF)])
Ginv=np.linalg.inv((omegaT+1j*ETA)*np.eye(nb)[None]-Hk_grid(kfg))
def gA(dm,dn): return np.einsum("k,kmn->mn",np.exp(2j*np.pi*(kfg[:,0]*dm+kfg[:,1]*dn)),Ginv)/NF**2
Gsub=np.zeros((dim,dim),complex)
for a,Ra in enumerate(sel):
    for b,Rb in enumerate(sel): Gsub[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=gA(Rg[Ra,0]-Rg[Rb,0],Rg[Ra,1]-Rg[Rb,1])
Tsub=np.linalg.solve(np.eye(dim)-Vsub@Gsub,Vsub).reshape(ns,nb,ns,nb)

# --- interpolate O^W to fine k_f (fixed k_i=K), band-project top valence ---
Hf=Hk_grid(kfg); ef,vf=np.linalg.eigh(Hf); topf=np.full(ef.shape[0],6)   # band 13 (VBM) at every k_f
vftop=np.take_along_axis(vf,topf[:,None,None],axis=2)[:,:,0]
# inverse FT consistent with toR (forward +k.R' on the bra R', -k.R on the ket R):
#   O^W(k_f,K) = sum_{R',R} e^{-2pi i k_f.R'} e^{+2pi i K.R} O^W(R',R)
phK=np.exp(2j*np.pi*(Rg@Kc)); Rsel=np.array([Rg[r] for r in sel]); phKsel=np.exp(2j*np.pi*(Rsel@Kc))
def bandmap(full=None,Tm=False):
    if Tm:
        A=np.einsum("aibj,b->aij",Tsub,phKsel); Okf=np.einsum("ka,aij->kij",np.exp(-2j*np.pi*(kfg@Rsel.T)),A)
    else:
        A=np.einsum("AiBj,B->Aij",full,phK);   Okf=np.einsum("kA,Aij->kij",np.exp(-2j*np.pi*(kfg@Rg.T)),A)
    return np.abs(np.einsum("ki,kij,j->k",vftop.conj(),Okf,vKtop))
mapM=bandmap(MWR); mapV=bandmap(VWR); mapT=bandmap(Tm=True)
iK=int(np.argmin(np.sum((kfg-Kc)**2,1)))
print(f"round-trip at k_f=K : |M|={mapM[iK]:.4f}  |V~|={mapV[iK]:.4f}  |T|={mapT[iK]:.4f} Ry")
print(f"  raw block band-13 : |M|=0.2464             |V~|=0.1121   (round-trip must match M,V~)")
print(f"max over BZ         : |M|={mapM.max():.4f}  |V~|={mapV.max():.4f}  |T|={mapT.max():.4f} Ry")

# --- Cartesian BZ ---
a1=np.array([240*0.150478,0.0])/6; a2=np.array([240*-0.075239,240*0.130318])/6
area=a1[0]*a2[1]-a1[1]*a2[0]; b1=2*np.pi/area*np.array([a2[1],-a2[0]]); b2=2*np.pi/area*np.array([-a1[1],a1[0]])
def fold(k1,k2):
    best=None
    for p in(-1,0,1):
        for q in(-1,0,1):
            c=(k1+p)*b1+(k2+q)*b2
            if best is None or c@c<best@best: best=c
    return best
kc=np.array([fold(kfg[k,0],kfg[k,1]) for k in range(len(kfg))]); Kpt=fold(1/3,1/3)
fig,ax=plt.subplots(1,3,figsize=(17,5.7))
panels=[(mapM,"$|M(k_f,K)|$  bare Born"),(mapV,r"$|\tilde V(k_f,K)|$  partial $T$"),(mapT,"$|T_{PP}(k_f,K)|$  full $T$")]
for p,(val,ttl) in enumerate(panels):
    sc=ax[p].tricontourf(kc[:,0],kc[:,1],val,levels=40,cmap="viridis")
    plt.colorbar(sc,ax=ax[p],label="Ry")
    ax[p].plot(*Kpt,"r*",ms=15,mew=1); ax[p].plot(0,0,"w+",ms=10,mew=2)
    ax[p].annotate("K",Kpt,xytext=(7,5),textcoords="offset points",color="r",fontweight="bold")
    ax[p].set_title(ttl+f"   (max {val.max():.3f} Ry)"); ax[p].set_xlabel("$k_x$ (1/Bohr)"); ax[p].set_ylabel("$k_y$"); ax[p].set_aspect("equal")
plt.suptitle(f"Wannier-interpolated 48×48 intraband maps, VBM (band 13), source K (★), on-shell ω=ε_VBM(K)={omegaT*RY:.2f} eV",y=1.01)
plt.tight_layout(); plt.savefig("p6_kmap_interp.png",dpi=130,bbox_inches="tight"); print("wrote p6_kmap_interp.png")
