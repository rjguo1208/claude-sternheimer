#!/usr/bin/env python3
"""On-shell diagonal T-matrix along the Gamma-M-K path (top valence band).
For each k on the path: omega = eps_top(k) (top valence band energy there); build the Koster-Slater
T^W = [1 - V~^W G^A(omega)]^{-1} V~^W with the fine-grid host G^A; Wannier-interpolate to (k_f=k_i=k)
and band-project the top band -> T_PP(k,k; eps_top(k)). Plot Re and Im vs the k-path."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05/RY; RCUT=4; NF=48; NSEG=90

fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]; eact=fb.read_reals(np.float64)
fb.read_record(np.complex128); fb.read_record(np.complex128)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5)); V4=V.reshape(nk,nb,nk,nb); Ek=eact.reshape(nk,nb)
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))

# V~^W(R',R) and H_W from the SAME filukk gauge
ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n) for m in ms for n in ms]); nR=len(Rg)
F=np.exp(2j*np.pi*(Rg@xk[:2]))
VW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),V4,Ub,optimize=True)
VWR=np.einsum("Af,fpgw,Bg->ApBw",F,VW,F.conj(),optimize=True)/nk**2
HWk=np.einsum("ksp,ks,ksw->kpw",Ub.conj(),Ek,Ub,optimize=True)
HWR=np.einsum("Rk,kpw->Rpw",np.exp(-2j*np.pi*(Rg@xk[:2])),HWk,optimize=True)/nk
def Hk1(kc2): return np.einsum("R,Rpw->pw",np.exp(2j*np.pi*(Rg@kc2)),HWR)

Rd=int(np.argmax([np.linalg.norm(VWR[r,:,r,:]) for r in range(nR)]))
def msh(R): return ((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2, ((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
sel=[r for r in range(nR) if max(abs(msh(r)[0]),abs(msh(r)[1]))<=RCUT]; ns=len(sel); dim=nb*ns
Rsel=np.array([Rg[r] for r in sel])
Vsub=np.zeros((dim,dim),complex)
for a in range(ns):
    for b in range(ns): Vsub[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=VWR[sel[a],:,sel[b],:]

# fine-grid eigendecomposition of H_W once; G^A(dR;omega) cheap per omega
kfg=np.array([[i/NF,j/NF] for i in range(NF) for j in range(NF)])
Hf=np.einsum("kR,Rpw->kpw",np.exp(2j*np.pi*(kfg@Rg.T)),HWR); efine,Wfine=np.linalg.eigh(Hf)
dRs=sorted({(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1]) for a in range(ns) for b in range(ns)})
dRidx={dR:i for i,dR in enumerate(dRs)}
phdR=np.array([np.exp(2j*np.pi*(kfg[:,0]*dm+kfg[:,1]*dn)) for (dm,dn) in dRs])
dRmat=np.array([[dRidx[(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1])] for b in range(ns)] for a in range(ns)])
def gA_table(omega):
    D=1.0/((omega+1j*ETA)-efine); G=np.einsum("kpn,kn,kqn->kpq",Wfine,D,Wfine.conj())
    return np.einsum("dk,kpq->dpq",phdR,G)/NF**2

def Tdiag(kc2):
    e,Vk=np.linalg.eigh(Hk1(kc2)); top=int(np.argmax(e)); om=e[top]; vtop=Vk[:,top]
    Gsub=gA_table(om)[dRmat].transpose(0,2,1,3).reshape(dim,dim)
    Tsub=np.linalg.solve(np.eye(dim)-Vsub@Gsub,Vsub).reshape(ns,nb,ns,nb)
    TWkk=np.einsum("a,aibj,b->ij",np.exp(-2j*np.pi*(Rsel@kc2)),Tsub,np.exp(2j*np.pi*(Rsel@kc2)))
    return vtop.conj()@TWkk@vtop, om

# path Gamma-M-K
Gpt=np.array([0.,0.]); Mpt=np.array([0.5,0.]); Kpt=np.array([1/3,1/3])
seg=lambda a,b,n:[a+(b-a)*t for t in np.linspace(0,1,n,endpoint=False)]
path=seg(Gpt,Mpt,NSEG)+seg(Mpt,Kpt,NSEG)+[Kpt]
a1=np.array([240*0.150478,0.0])/6; a2=np.array([240*-0.075239,240*0.130318])/6
area=a1[0]*a2[1]-a1[1]*a2[0]; b1=2*np.pi/area*np.array([a2[1],-a2[0]]); b2=2*np.pi/area*np.array([-a1[1],a1[0]])
kc=np.array([k[0]*b1+k[1]*b2 for k in path])
dist=np.concatenate([[0],np.cumsum(np.linalg.norm(np.diff(kc,axis=0),axis=1))])
xM,xK=dist[NSEG],dist[-1]
Re=np.empty(len(path)); Im=np.empty(len(path)); eps=np.empty(len(path))
for i,k in enumerate(path):
    t,om=Tdiag(np.array(k)); Re[i],Im[i],eps[i]=t.real,t.imag,om
print(f"Gamma: eps={eps[0]*RY:.3f} eV  T={Re[0]:+.4f}{Im[0]:+.4f}i")
print(f"M    : eps={eps[NSEG]*RY:.3f} eV  T={Re[NSEG]:+.4f}{Im[NSEG]:+.4f}i")
print(f"K(VBM): eps={eps[-1]*RY:.3f} eV  T={Re[-1]:+.4f}{Im[-1]:+.4f}i  |T|={abs(Re[-1]+1j*Im[-1]):.4f}")

fig,ax=plt.subplots(figsize=(8.4,5.0))
ax.axhline(0,c="0.6",lw=0.6)
ax.plot(dist,Re,c="#1f77b4",lw=1.9,label=r"Re $T_{PP}(k,k;\varepsilon_k)$")
ax.plot(dist,Im,c="#d62728",lw=1.9,label=r"Im $T_{PP}(k,k;\varepsilon_k)$")
ax.plot(dist,np.abs(Re+1j*Im),c="k",lw=1.0,ls=":",label=r"$|T_{PP}|$")
for x in (0,xM,xK): ax.axvline(x,c="0.8",lw=0.8)
ax.set_xticks([0,xM,xK]); ax.set_xticklabels(["Γ","M","K"]); ax.set_xlim(0,xK)
ax.set_xlabel("k-path"); ax.set_ylabel(r"$T_{PP}(k,k;\varepsilon_{\rm top}(k))$  (Ry)")
ax.set_title("On-shell diagonal $T$-matrix along Γ–M–K (top valence band)")
ax.legend(loc="upper left"); ax.grid(alpha=.3)
axb=ax.twinx(); axb.plot(dist,eps*RY,c="0.55",lw=1.2,ls="--"); axb.set_ylabel(r"$\varepsilon_{\rm top}(k)$  (eV)",color="0.45")
axb.tick_params(axis="y",colors="0.45")
plt.tight_layout(); plt.savefig("p6_tpath.png",dpi=130); print("wrote p6_tpath.png")
