#!/usr/bin/env python3
"""On-shell diagonal T-matrix along Gamma-M-K (top valence band), WITH vs WITHOUT rest-space.
For each k: omega = eps_top(k); Koster-Slater T_X = [1 - X^W G^A(omega)]^{-1} X^W with the fine-grid
host G^A, for X = V~ (with rest-space) and X = M (no rest-space); Wannier-interpolate to (k_f=k_i=k)
and band-project the top band. Plot Re and Im of both vs the k-path."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05/RY; RCUT=4; NF=48; NSEG=90

fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]; eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")     # bare matrix element (no rest-space)
fb.read_record(np.complex128)                                    # Sigma (skip)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()   # V~ = M + Sigma
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
MWR=toR(M4); VWR=toR(V4)
HWk=np.einsum("ksp,ks,ksw->kpw",Ub.conj(),Ek,Ub,optimize=True)
HWR=np.einsum("Rk,kpw->Rpw",np.exp(-2j*np.pi*(Rg@xk[:2])),HWk,optimize=True)/nk
def Hk1(kc2): return np.einsum("R,Rpw->pw",np.exp(2j*np.pi*(Rg@kc2)),HWR)

Rd=int(np.argmax([np.linalg.norm(VWR[r,:,r,:]) for r in range(nR)]))
def msh(R): return ((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2, ((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
sel=[r for r in range(nR) if max(abs(msh(r)[0]),abs(msh(r)[1]))<=RCUT]; ns=len(sel); dim=nb*ns
Rsel=np.array([Rg[r] for r in sel])
def sub(OWR):
    S=np.zeros((dim,dim),complex)
    for a in range(ns):
        for b in range(ns): S[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=OWR[sel[a],:,sel[b],:]
    return S
Vsub=sub(VWR); Msub=sub(MWR)

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
    Gsub=gA_table(om)[dRmat].transpose(0,2,1,3).reshape(dim,dim); I=np.eye(dim)
    phb=np.exp(-2j*np.pi*(Rsel@kc2)); phk=np.exp(2j*np.pi*(Rsel@kc2))
    def proj(Xsub):
        Tsub=np.linalg.solve(I-Xsub@Gsub,Xsub).reshape(ns,nb,ns,nb)
        return vtop.conj()@np.einsum("a,aibj,b->ij",phb,Tsub,phk)@vtop
    return proj(Vsub), proj(Msub), om

Gpt=np.array([0.,0.]); Mpt=np.array([0.5,0.]); Kpt=np.array([1/3,1/3])
seg=lambda a,b,n:[a+(b-a)*t for t in np.linspace(0,1,n,endpoint=False)]
path=seg(Gpt,Mpt,NSEG)+seg(Mpt,Kpt,NSEG)+[Kpt]
a1=np.array([240*0.150478,0.0])/6; a2=np.array([240*-0.075239,240*0.130318])/6
area=a1[0]*a2[1]-a1[1]*a2[0]; b1=2*np.pi/area*np.array([a2[1],-a2[0]]); b2=2*np.pi/area*np.array([-a1[1],a1[0]])
kc=np.array([k[0]*b1+k[1]*b2 for k in path]); dist=np.concatenate([[0],np.cumsum(np.linalg.norm(np.diff(kc,axis=0),axis=1))])
xM,xK=dist[NSEG],dist[-1]
Tv=np.empty(len(path),complex); Tm=np.empty(len(path),complex); eps=np.empty(len(path))
for i,k in enumerate(path):
    tv,tm,om=Tdiag(np.array(k)); Tv[i]=tv; Tm[i]=tm; eps[i]=om
    if i%30==0: print(f"  {i}/{len(path)} done",flush=True)
for nm,j in [("Gamma",0),("M",NSEG),("K(VBM)",-1)]:
    print(f"{nm}: eps={eps[j]*RY:.3f} eV  T_PP(with rest)={Tv[j]:+.4f}  T_M(no rest)={Tm[j]:+.4f}")

fig,ax=plt.subplots(figsize=(8.8,5.3))
ax.axhline(0,c="0.6",lw=0.6)
ax.plot(dist,Tv.real,c="#1f77b4",lw=1.9,label=r"Re $T_{PP}$ (with rest-space)")
ax.plot(dist,Tv.imag,c="#d62728",lw=1.9,label=r"Im $T_{PP}$ (with rest-space)")
ax.plot(dist,Tm.real,c="#1f77b4",lw=1.5,ls="--",label=r"Re $T_M$ (no rest-space)")
ax.plot(dist,Tm.imag,c="#d62728",lw=1.5,ls="--",label=r"Im $T_M$ (no rest-space)")
for x in (0,xM,xK): ax.axvline(x,c="0.8",lw=0.8)
ax.set_xticks([0,xM,xK]); ax.set_xticklabels(["Γ","M","K"]); ax.set_xlim(0,xK)
ax.set_xlabel("k-path"); ax.set_ylabel(r"$T(k,k;\varepsilon_{\rm top}(k))$  (Ry)")
ax.set_title("On-shell diagonal $T$-matrix along Γ–M–K: with vs without rest-space")
ax.legend(loc="upper left",fontsize=8,ncol=2); ax.grid(alpha=.3)
axb=ax.twinx(); axb.plot(dist,eps*RY,c="0.6",lw=1.1,ls=":"); axb.set_ylabel(r"$\varepsilon_{\rm top}(k)$  (eV)",color="0.45"); axb.tick_params(axis="y",colors="0.45")
plt.tight_layout(); plt.savefig("p6_tpath.png",dpi=130); print("wrote p6_tpath.png")
