#!/usr/bin/env python3
"""Electron-defect self-energy and spectral function along Gamma-M-K (VBM, band 13).
  Sigma_e-def(nk,omega) = n_d * T_{nk,nk}(omega)            (single-site / dilute, full T with rest-space)
  Dyson:  G(nk,omega) = 1/(omega - eps_nk - Sigma(nk,omega))
  Spectral fn:  A(nk,omega) = -(1/pi) Im G = (1/pi)(-Im Sigma)/[(omega-eps_nk-Re Sigma)^2 + (Im Sigma)^2]
n_d = defects per primitive cell (argv[1], default 0.01 ~ 1.2e13 cm^-2 for MoS2).
Koster-Slater T_sub(omega) is k-independent -> one solve per omega; band 13 = eigh sorted slot 6."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05/RY; RCUT=4; NF=48; NSEG=70; NW=260; VBM=6
ND=float(sys.argv[1]) if len(sys.argv)>1 else 0.01          # defects per primitive cell
INPUT=sys.argv[2] if len(sys.argv)>2 else "V"              # V=with rest-space (V~=M+Sigma), M=no rest-space (bare M)
restlbl="with rest-space" if INPUT=="V" else "no rest-space (bare $M$)"
ACELL=8.64e-16                                              # MoS2 primitive cell area (cm^2)

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
VWR=toR(V4); XWR=toR(M4) if INPUT=="M" else VWR      # Rd/sel from V~; Koster-Slater input = X
HWk=np.einsum("ksp,ks,ksw->kpw",Ub.conj(),Ek,Ub,optimize=True)
HWR=np.einsum("Rk,kpw->Rpw",np.exp(-2j*np.pi*(Rg@xk[:2])),HWk,optimize=True)/nk
def Hk1(k): return np.einsum("R,Rpw->pw",np.exp(2j*np.pi*(Rg@k)),HWR)

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
    e,Vk=np.linalg.eigh(Hk1(k)); eps[i]=e[VBM]; vt=Vk[:,VBM]
    phb=np.exp(-2j*np.pi*(Rsel@k)); phk=np.exp(2j*np.pi*(Rsel@k))
    lvec[i]=(phb[:,None]*vt.conj()[None,:]).reshape(dim); rvec[i]=(phk[:,None]*vt[None,:]).reshape(dim)
a1=np.array([240*0.150478,0.0])/6; a2=np.array([240*-0.075239,240*0.130318])/6
area=a1[0]*a2[1]-a1[1]*a2[0]; b1=2*np.pi/area*np.array([a2[1],-a2[0]]); b2=2*np.pi/area*np.array([-a1[1],a1[0]])
kc=np.array([k[0]*b1+k[1]*b2 for k in path]); dist=np.concatenate([[0],np.cumsum(np.linalg.norm(np.diff(kc,axis=0),axis=1))])
xM,xK=dist[NSEG],dist[-1]

WLO=eps.min()*RY-0.8; WHI=eps.max()*RY+0.9
omega=np.linspace(WLO,WHI,NW)/RY; T=np.empty((Nk,NW),complex)
for j,om in enumerate(omega):
    T[:,j]=((lvec@Tsub_of(om))*rvec).sum(1)
    if j%40==0: print(f"  {j}/{NW}",flush=True)

# self-energy and Dyson spectral function (work in eV)
Sig=ND*T*RY                                  # Sigma(nk,omega) in eV
om_eV=omega*RY; eps_eV=eps*RY
denom=(om_eV[None,:]-eps_eV[:,None]-Sig.real)**2+Sig.imag**2
A=(1/np.pi)*(-Sig.imag)/denom                # spectral function (1/eV), >=0
print(f"n_d={ND} ({ND/ACELL:.1e} cm^-2);  max -Im Sigma={-Sig.imag.min()*1e3:.1f} meV;  max A={A.max():.1f} /eV")

fig,ax=plt.subplots(1,3,figsize=(16.5,5.2))
# Re Sigma (meV) -- signed -> symlog
ReS=Sig.real.T*1e3; vR=np.percentile(np.abs(ReS),99.5)
p0=ax[0].pcolormesh(dist,om_eV,ReS,shading="gouraud",cmap="RdBu_r",
                    norm=SymLogNorm(linthresh=max(vR/40,1e-2),vmin=-vR,vmax=vR,base=10))
plt.colorbar(p0,ax=ax[0],label="meV"); ax[0].set_title(r"Re $\Sigma(nk,\omega)$  (level shift, symlog)")
# -Im Sigma (meV) >=0 -> log
mIm=-Sig.imag.T*1e3; vI=np.percentile(mIm,99.7); mIm=np.clip(mIm,vI/1e2,None)
p1=ax[1].pcolormesh(dist,om_eV,mIm,shading="gouraud",cmap="magma",norm=LogNorm(vmin=vI/1e2,vmax=vI))
plt.colorbar(p1,ax=ax[1],label="meV"); ax[1].set_title(r"$-$Im $\Sigma(nk,\omega)$  (rate $\Gamma/2$, log)")
# A(nk,omega) >=0 -> log (spans the QP peak down to the faint in-gap/tail weight)
AA=A.T; vA=np.percentile(AA,99.8); AA=np.clip(AA,vA/1e3,None)
p2=ax[2].pcolormesh(dist,om_eV,AA,shading="gouraud",cmap="magma",norm=LogNorm(vmin=vA/1e3,vmax=vA))
plt.colorbar(p2,ax=ax[2],label="1/eV"); ax[2].set_title(r"$A(nk,\omega)$  spectral function (log)")
for p in (0,1,2):
    ax[p].plot(dist,eps_eV,c=("k" if p==0 else "w"),lw=1.0,ls="--")   # bare band eps_nk
    ax[p].axvline(xM,c=("0.5" if p==0 else "w"),lw=0.6,alpha=0.6)
    ax[p].set_xticks([0,xM,xK]); ax[p].set_xticklabels(["Γ","M","K"]); ax[p].set_xlim(0,xK)
    ax[p].set_ylim(WLO,WHI); ax[p].set_xlabel("k-path"); ax[p].set_ylabel(r"$\omega$  (eV)")
plt.suptitle(f"e-defect $\\Sigma=n_d T$ and Dyson spectral function $A$ on the VBM (band 13), {restlbl}, $n_d={ND}$ ({ND/ACELL:.0e} cm$^{{-2}}$);  LOG colour scale;  dashed = bare $\\varepsilon_{{nk}}$",y=1.02)
out=f"p6_spectral_{INPUT}.png"; plt.tight_layout(); plt.savefig(out,dpi=130,bbox_inches="tight"); print("wrote",out)
