#!/usr/bin/env python3
"""Diagonal T-matrix T(nk,omega) color map along Gamma-M-K vs energy (top valence band).
Self-energy Sigma_e-def(nk,omega) = n_d * T_{nk,nk}(omega).
  argv[1] = "V" (default) -> full T with rest-space (input V~ = M + Sigma)
          = "M"           -> no rest-space (bare matrix element M)
Key: the Koster-Slater T_sub(omega)=[1-X^W G^A(omega)]^{-1}X^W is k-INDEPENDENT (only X and the host
G^A(omega) enter) -> solve once per omega, then T(nk,omega)=l(k).T_sub(omega).r(k) for all path-k via
the Wannier interpolation + top-band projection. Re (diverging) and -Im (sequential) maps."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05/RY; RCUT=4; NF=48; NSEG=70; NW=141; WLO,WHI=-3.5,-0.5
INPUT=sys.argv[1] if len(sys.argv)>1 else "V"
tag="no rest-space (bare $M$)" if INPUT=="M" else "full $T$ with rest-space"

fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]; eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.read_record(np.complex128)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5)); M4=M.reshape(nk,nb,nk,nb); V4=V.reshape(nk,nb,nk,nb); Ek=eact.reshape(nk,nb)
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))

ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n) for m in ms for n in ms]); nR=len(Rg)
F=np.exp(2j*np.pi*(Rg@xk[:2]))
def toR(O4):
    OW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),O4,Ub,optimize=True)
    return np.einsum("Af,fpgw,Bg->ApBw",F,OW,F.conj(),optimize=True)/nk**2
VWR=toR(V4); XWR=toR(M4) if INPUT=="M" else VWR     # Rd/sel from V~; Koster-Slater input = X
HWk=np.einsum("ksp,ks,ksw->kpw",Ub.conj(),Ek,Ub,optimize=True)
HWR=np.einsum("Rk,kpw->Rpw",np.exp(-2j*np.pi*(Rg@xk[:2])),HWk,optimize=True)/nk
def Hk1(kc2): return np.einsum("R,Rpw->pw",np.exp(2j*np.pi*(Rg@kc2)),HWR)

Rd=int(np.argmax([np.linalg.norm(VWR[r,:,r,:]) for r in range(nR)]))
def msh(R): return ((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2, ((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
sel=[r for r in range(nR) if max(abs(msh(r)[0]),abs(msh(r)[1]))<=RCUT]; ns=len(sel); dim=nb*ns
Rsel=np.array([Rg[r] for r in sel])
Xsub=np.zeros((dim,dim),complex)
for a in range(ns):
    for b in range(ns): Xsub[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=XWR[sel[a],:,sel[b],:]

kfg=np.array([[i/NF,j/NF] for i in range(NF) for j in range(NF)])
Hf=np.einsum("kR,Rpw->kpw",np.exp(2j*np.pi*(kfg@Rg.T)),HWR); efine,Wfine=np.linalg.eigh(Hf)
dRs=sorted({(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1]) for a in range(ns) for b in range(ns)})
dRidx={dR:i for i,dR in enumerate(dRs)}
phdR=np.array([np.exp(2j*np.pi*(kfg[:,0]*dm+kfg[:,1]*dn)) for (dm,dn) in dRs])
dRmat=np.array([[dRidx[(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1])] for b in range(ns)] for a in range(ns)])
def Tsub_of(omega):
    D=1.0/((omega+1j*ETA)-efine); G=np.einsum("kpn,kn,kqn->kpq",Wfine,D,Wfine.conj())
    gt=(np.einsum("dk,kpq->dpq",phdR,G)/NF**2)[dRmat].transpose(0,2,1,3).reshape(dim,dim)
    return np.linalg.solve(np.eye(dim)-Xsub@gt,Xsub)

Gpt=np.array([0.,0.]); Mpt=np.array([0.5,0.]); Kpt=np.array([1/3,1/3])
seg=lambda a,b,n:[a+(b-a)*t for t in np.linspace(0,1,n,endpoint=False)]
path=np.array(seg(Gpt,Mpt,NSEG)+seg(Mpt,Kpt,NSEG)+[Kpt]); Nk=len(path)
lvec=np.empty((Nk,dim),complex); rvec=np.empty((Nk,dim),complex); eps=np.empty(Nk)
for i,k in enumerate(path):
    e,Vk=np.linalg.eigh(Hk1(k)); top=int(np.argmax(e)); eps[i]=e[top]; vt=Vk[:,top]
    phb=np.exp(-2j*np.pi*(Rsel@k)); phk=np.exp(2j*np.pi*(Rsel@k))
    lvec[i]=(phb[:,None]*vt.conj()[None,:]).reshape(dim); rvec[i]=(phk[:,None]*vt[None,:]).reshape(dim)
a1=np.array([240*0.150478,0.0])/6; a2=np.array([240*-0.075239,240*0.130318])/6
area=a1[0]*a2[1]-a1[1]*a2[0]; b1=2*np.pi/area*np.array([a2[1],-a2[0]]); b2=2*np.pi/area*np.array([-a1[1],a1[0]])
kc=np.array([k[0]*b1+k[1]*b2 for k in path]); dist=np.concatenate([[0],np.cumsum(np.linalg.norm(np.diff(kc,axis=0),axis=1))])
xM,xK=dist[NSEG],dist[-1]

omega=np.linspace(WLO,WHI,NW)/RY; Tmap=np.empty((Nk,NW),complex)
for j,om in enumerate(omega):
    Tmap[:,j]=((lvec@Tsub_of(om))*rvec).sum(1)
    if j%30==0: print(f"  {j}/{NW}",flush=True)
iK=Nk-1; jVBM=int(np.argmin(np.abs(omega-omega0)))
ref="-0.033-0.081i (T_M)" if INPUT=="M" else "-0.033-0.097i (T_PP)"
print(f"[{INPUT}] check T(K, omega~VBM)={Tmap[iK,jVBM]:+.4f}  (Figure-3 on-shell {ref})")
print(f"Re range [{Tmap.real.min():.3f},{Tmap.real.max():.3f}]  Im range [{Tmap.imag.min():.3f},{Tmap.imag.max():.3f}]")

fig,ax=plt.subplots(1,2,figsize=(14,5.4))
vR=np.percentile(np.abs(Tmap.real),99)
p0=ax[0].pcolormesh(dist,omega*RY,Tmap.real.T,shading="gouraud",cmap="RdBu_r",vmin=-vR,vmax=vR)
plt.colorbar(p0,ax=ax[0],label="Ry"); ax[0].plot(dist,eps*RY,c="k",lw=1.4,ls="--")
ax[0].set_title(r"Re $T(nk,\omega)$  (level shift)")
mIm=-Tmap.imag; vI=np.percentile(mIm,99)
p1=ax[1].pcolormesh(dist,omega*RY,mIm.T,shading="gouraud",cmap="magma",vmin=0,vmax=vI)
plt.colorbar(p1,ax=ax[1],label="Ry"); ax[1].plot(dist,eps*RY,c="w",lw=1.4,ls="--")
ax[1].set_title(r"$-$Im $T(nk,\omega)$  (spectral weight $\propto$ rate)")
for p in (0,1):
    ax[p].axvline(xM,c="0.5" if p==0 else "w",lw=0.7,alpha=0.6)
    ax[p].set_xticks([0,xM,xK]); ax[p].set_xticklabels(["Γ","M","K"]); ax[p].set_xlim(0,xK)
    ax[p].set_ylim(WLO,WHI); ax[p].set_xlabel("k-path"); ax[p].set_ylabel(r"$\omega$  (eV)")
plt.suptitle(r"Diagonal $T(nk,\omega)$ (top valence band, "+tag+r") along Γ–M–K — dashed = on-shell $\varepsilon_{\rm top}(k)$",y=1.02)
out="p6_tmap_M.png" if INPUT=="M" else "p6_tmap.png"
plt.tight_layout(); plt.savefig(out,dpi=130,bbox_inches="tight"); print("wrote",out)
