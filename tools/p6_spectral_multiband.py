#!/usr/bin/env python3
"""Complete MULTIBAND electron-defect spectral function along Gamma-M-K (all active bands).
  Sigma_{nn'}(k,omega) = n_d T_{(nk),(n'k)}(omega)          (band-space matrix from the full T)
  G(k,omega) = [ omega I - H0(k) - Sigma(k,omega) ]^{-1},   H0=diag(eps_nk)
  A(k,omega) = -(1/pi) Im Tr G                              (the trace sums all bands)
Computed for with-rest (V~) and no-rest (M), n_d = 1% and 5%, log colour scale; SAME energy window
as the single-band figures (band-13 auto window).  T_sub(omega) is k- and band-independent (one solve
per omega per input); we keep the full Wannier T^W(k,k) and rotate to all bands V(k)^dag T^W V(k)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05/RY; RCUT=4; NF=48; NSEG=50; NW=280; VBM=6; CBM=7; ACELL=8.64e-16
NDS=[0.01,0.05]
import sys
BLOCK=sys.argv[1] if len(sys.argv)>1 else "vtilde_block.dat"   # which Vtilde block (reference)
TAG  =sys.argv[2] if len(sys.argv)>2 else "wide"               # output-name suffix

fb=FortranFile(BLOCK,"r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]; eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.read_record(np.complex128)
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
VWR=toR(V4); MWR=toR(M4)
HWk=np.einsum("ksp,ks,ksw->kpw",Ub.conj(),Ek,Ub,optimize=True)
HWR=np.einsum("Rk,kpw->Rpw",np.exp(-2j*np.pi*(Rg@xk[:2])),HWk,optimize=True)/nk
def Hk1(k): return np.einsum("R,Rpw->pw",np.exp(2j*np.pi*(Rg@k)),HWR)

Rd=int(np.argmax([np.linalg.norm(VWR[r,:,r,:]) for r in range(nR)]))
def msh(R): return ((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2, ((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
sel=[r for r in range(nR) if max(abs(msh(r)[0]),abs(msh(r)[1]))<=RCUT]; ns=len(sel); dim=nb*ns
Rsel=np.array([Rg[r] for r in sel])
def sub(WR):
    S=np.zeros((dim,dim),complex)
    for a in range(ns):
        for b in range(ns): S[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=WR[sel[a],:,sel[b],:]
    return S
Vsub=sub(VWR); Msub=sub(MWR)
kfg=np.array([[i/NF,j/NF] for i in range(NF) for j in range(NF)])
Hf=np.einsum("kR,Rpw->kpw",np.exp(2j*np.pi*(kfg@Rg.T)),HWR); efine,Wfine=np.linalg.eigh(Hf)
dRs=sorted({(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1]) for a in range(ns) for b in range(ns)})
dRidx={dR:i for i,dR in enumerate(dRs)}
phdR=np.array([np.exp(2j*np.pi*(kfg[:,0]*dm+kfg[:,1]*dn)) for (dm,dn) in dRs])
dRmat=np.array([[dRidx[(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1])] for b in range(ns)] for a in range(ns)])
def Tsub_of(omega,Xsub):
    D=1.0/((omega+1j*ETA)-efine); G=np.einsum("kpn,kn,kqn->kpq",Wfine,D,Wfine.conj())
    gt=(np.einsum("dk,kpq->dpq",phdR,G)/NF**2)[dRmat].transpose(0,2,1,3).reshape(dim,dim)
    return np.linalg.solve(np.eye(dim)-Xsub@gt,Xsub)

Gpt=np.array([0.,0.]); Mpt=np.array([0.5,0.]); Kpt=np.array([1/3,1/3])
seg=lambda a,b,n:[a+(b-a)*t for t in np.linspace(0,1,n,endpoint=False)]
path=np.array(seg(Gpt,Mpt,NSEG)+seg(Mpt,Kpt,NSEG)+[Kpt]); Nk=len(path)
Vall=np.empty((Nk,nb,nb),complex); EPS=np.empty((Nk,nb)); phbL=np.empty((Nk,ns),complex); phkL=np.empty((Nk,ns),complex)
for i,k in enumerate(path):
    e,Vk=np.linalg.eigh(Hk1(k)); EPS[i]=e; Vall[i]=Vk
    phbL[i]=np.exp(-2j*np.pi*(Rsel@k)); phkL[i]=np.exp(2j*np.pi*(Rsel@k))
a1=np.array([240*0.150478,0.0])/6; a2=np.array([240*-0.075239,240*0.130318])/6
area=a1[0]*a2[1]-a1[1]*a2[0]; b1=2*np.pi/area*np.array([a2[1],-a2[0]]); b2=2*np.pi/area*np.array([-a1[1],a1[0]])
kc=np.array([k[0]*b1+k[1]*b2 for k in path]); dist=np.concatenate([[0],np.cumsum(np.linalg.norm(np.diff(kc,axis=0),axis=1))])
xM,xK=dist[NSEG],dist[-1]
WLO=EPS[:,VBM].min()*RY-0.8; WHI=-3.0       # lower edge unchanged; upper extended to -3.0 eV (into the conduction manifold)
omega=np.linspace(WLO,WHI,NW)/RY; om_eV=omega*RY

def Tband(Xsub,tag):                                          # [Nk,NW,nb,nb]  band-space T_{nn'}(k;omega)
    Tb=np.empty((Nk,NW,nb,nb),complex)
    for j,om in enumerate(omega):
        Ts=Tsub_of(om,Xsub).reshape(ns,nb,ns,nb)
        TW=np.einsum("ka,apbq,kb->kpq",phbL,Ts,phkL,optimize=True)            # T^W(k,k) [Nk,nb,nb]
        Tb[:,j]=np.einsum("kpn,kpq,kqm->knm",Vall.conj(),TW,Vall,optimize=True)
        if j%40==0: print(f"  {tag} {j}/{NW}",flush=True)
    return Tb
TbV=Tband(Vsub,"V"); TbM=Tband(Msub,"M")

Inb=np.eye(nb); H0d=(EPS*RY)[:,:,None]*Inb[None]              # diag(eps_nk) eV  [Nk,nb,nb]
def Aspec(Tb,nd):
    Sig=nd*Tb*RY
    Am=om_eV[None,:,None,None]*Inb[None,None]-H0d[:,None]-Sig
    return -(1/np.pi)*np.imag(np.einsum("kjnn->kj",np.linalg.inv(Am)))        # A(k,omega) 1/eV
maps={(nd,inp):Aspec(Tb,nd) for nd in NDS for inp,Tb in (("V",TbV),("M",TbM))}
np.savez(f"p6_multiband_maps_{TAG}.npz",dist=dist,om_eV=om_eV,EPS=EPS,xM=xM,xK=xK,WLO=WLO,WHI=WHI,
         **{f"{int(nd*100)}_{inp}":maps[(nd,inp)] for nd in NDS for inp in ("V","M")})
allA=np.concatenate([m.ravel() for m in maps.values()]); vmax=np.percentile(allA,99.8); vmin=vmax/1e3
gap=(EPS[:,CBM].min()-EPS[:,VBM].max())*RY
print(f"window [{WLO:.2f},{WHI:.2f}] eV; VBM={EPS[:,VBM].max()*RY:.3f}, CBM={EPS[:,CBM].min()*RY:.3f}, gap={gap:.3f} eV")
print(f"bands in window at K: {[7+s for s in range(nb) if WLO<=EPS[-1,s]*RY<=WHI]};  A range [{allA.min():.2g},{allA.max():.1f}] /eV")

fig,ax=plt.subplots(2,2,figsize=(13,9))
for r,nd in enumerate(NDS):
    for c,(inp,lbl) in enumerate([("V","with rest-space ($\\tilde V$)"),("M","no rest-space ($M$)")]):
        A=maps[(nd,inp)]; a=ax[r,c]
        pc=a.pcolormesh(dist,om_eV,np.clip(A.T,vmin,None),shading="gouraud",cmap="magma",norm=LogNorm(vmin=vmin,vmax=vmax))
        for s in range(nb): a.plot(dist,EPS[:,s]*RY,c="c",lw=0.6,ls="--",alpha=0.55)   # bare bands
        a.axvline(xM,c="w",lw=0.6,alpha=0.5); a.set_xticks([0,xM,xK]); a.set_xticklabels(["Γ","M","K"])
        a.set_xlim(0,xK); a.set_ylim(WLO,WHI); a.set_ylabel(r"$\omega$ (eV)")
        a.set_title(f"$n_d={int(nd*100)}\\%$, {lbl}",fontsize=10)
        if r==1: a.set_xlabel("k-path")
        plt.colorbar(pc,ax=a,label="A  (1/eV)")
plt.suptitle(r"Complete multiband spectral function  $A(k,\omega)=-\frac{1}{\pi}\,\mathrm{Im}\,\mathrm{Tr}\,[\omega-H_0-n_d\,T]^{-1}$  (all active bands; cyan = bare $\varepsilon_{nk}$; log scale)",y=1.0)
plt.tight_layout(); plt.savefig(f"p6_spectral_multiband_{TAG}.png",dpi=130,bbox_inches="tight"); print(f"wrote p6_spectral_multiband_{TAG}.png  (window [{WLO:.2f},{WHI:.2f}] eV, ref omega0={omega0*RY:.3f} eV)")
