#!/usr/bin/env python3
"""On-shell defect-scattering rates for states near the VBM and CBM, comparing
  (a) 11-band active space, bare M        (no rest dressing)
  (b) 11-band + FULL-ORDER rest Sigma     (fesh60 direct-resolvent Vtilde)
Rate:  hbar/tau_nk = -2 n_d Im T_nn(k, eps_nk + i eta)   (optical theorem, on-shell diag T)
T from the Wannier real-space sub-block (R_cut=4) + FINE-grid host GF (NF=48), exactly the
machinery of p6_spectral_multiband.py; both M and Vtilde records come from the SAME block file.
States = all fine-grid (n,k) with eps within WIN of the band edge. n_d = 1% (linear scaling)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05/RY; RCUT=4; NF=48; WIN=0.30; ND=0.01; HBAR=658.2119569  # meV*fs

fb=FortranFile("vtilde_block_fesh60.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
fb.read_ints(np.int32); eact=fb.read_reals(np.float64)
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
def Tsub_of(omega,Xsub,born=False):
    D=1.0/((omega+1j*ETA)-efine); G=np.einsum("kpn,kn,kqn->kpq",Wfine,D,Wfine.conj())
    gt=(np.einsum("dk,kpq->dpq",phdR,G)/NF**2)[dRmat].transpose(0,2,1,3).reshape(dim,dim)
    if born: return Xsub+Xsub@gt@Xsub                 # 1st Born: golden rule |X|^2 via optical theorem
    return np.linalg.solve(np.eye(dim)-Xsub@gt,Xsub)

e_eV=efine*RY                                       # [NF^2, nb] fine-grid band energies (eV)
VBMe=e_eV[:,6].max(); CBMe=e_eV[:,7].min()
print(f"fine-grid edges: VBM={VBMe:.4f}  CBM={CBMe:.4f}  gap={CBMe-VBMe:.3f} eV",flush=True)
states=[]                                           # (ik, n, eps, side)
for n in range(nb):
    for ik in range(NF*NF):
        e=e_eV[ik,n]
        if VBMe-WIN<=e<=VBMe+1e-9: states.append((ik,n,e,"V"))
        elif CBMe-1e-9<=e<=CBMe+WIN: states.append((ik,n,e,"C"))
print(f"selected {len(states)} states (|eps-edge|<={WIN} eV); bands: "
      f"{sorted({n for _,n,_,_ in states})}",flush=True)
ksel=sorted({ik for ik,_,_,_ in states}); kpos={ik:i for i,ik in enumerate(ksel)}
phb=np.exp(-2j*np.pi*(kfg[ksel]@Rsel[:,:2].T)); phk=phb.conj()       # [Nsel, ns]
Wsel=Wfine[ksel]                                                     # [Nsel, nb, nb]

NWG=40
omV=np.linspace(VBMe-WIN-0.02,VBMe+0.02,NWG); omC=np.linspace(CBMe-0.02,CBMe+WIN+0.02,NWG)
omall=np.concatenate([omV,omC])
def diagT(Xsub,tag,born=False):
    out=np.empty((len(omall),len(ksel),nb))          # Im T_nn(k;omega)  (Ry)
    for j,om in enumerate(omall):
        Ts=Tsub_of(om/RY,Xsub,born).reshape(ns,nb,ns,nb)
        TW=np.einsum("ka,apbq,kb->kpq",phb,Ts,phk,optimize=True)
        out[j]=np.imag(np.einsum("kpn,kpq,kqn->kn",Wsel.conj(),TW,Wsel,optimize=True))
        if j%10==0: print(f"  {tag} {j}/{len(omall)}",flush=True)
    return out
ImM=diagT(Msub,"M  (bare)"); ImV=diagT(Vsub,"V~ (full)"); ImB=diagT(Msub,"|M|^2 Born",born=True)

def rates(Im):                                       # hbar/tau in meV at n_d=ND
    r=np.empty(len(states))
    for i,(ik,n,e,side) in enumerate(states):
        og=omV if side=="V" else omC; sl=slice(0,NWG) if side=="V" else slice(NWG,2*NWG)
        r[i]=-2.0*ND*np.interp(e,og,Im[sl,kpos[ik],n])*RY*1e3
    return r
rM=rates(ImM); rV=rates(ImV); rB=rates(ImB)
eps=np.array([s[2] for s in states]); side=np.array([s[3] for s in states])
np.savez("p7_scatrate.npz",eps=eps,side=side,rM=rM,rV=rV,rB=rB,VBM=VBMe,CBM=CBMe,nd=ND,eta=ETA*RY,
         band=np.array([s[1] for s in states]),ik=np.array([s[0] for s in states]))

fig,axs=plt.subplots(1,2,figsize=(12.4,5.2),sharey=True)
for ax,sd,edge,lbl in [(axs[0],"V",VBMe,"VBM"),(axs[1],"C",CBMe,"CBM")]:
    m=side==sd; x=eps[m]-edge
    ax.scatter(x,rB[m],s=10,c="#d17b0f",alpha=.6,label=r"golden rule, bare $|M|^2$ (1st Born)")
    ax.scatter(x,rM[m],s=10,c="#2a7a2a",alpha=.6,label="full $T$, bare $M$ (11 bands)")
    ax.scatter(x,rV[m],s=10,c="#6a4c93",alpha=.6,label=r"full $T$, full-order $\tilde V$ (fesh60)")
    ax.set_yscale("log"); ax.axvline(0,color="#999",lw=.8,ls=":")
    ax.set_xlabel(rf"$\varepsilon_{{n\mathbf{{k}}}}-${lbl} (eV)")
    ax.set_title(f"states near the {lbl}",fontsize=11)
    ax.grid(alpha=.25,which="both")
axs[0].set_ylabel(r"$\hbar/\tau_{n\mathbf{k}}$ (meV)  at  $n_d=1\%$")
axs[0].invert_xaxis()
h,l=axs[0].get_legend_handles_labels()
fig.legend(h,l,ncol=3,loc="lower center",bbox_to_anchor=(0.5,0.90),fontsize=9.5,frameon=False)
fig.suptitle(r"Defect-limited scattering rate  $\hbar/\tau_{n\mathbf{k}}=-2 n_d\,\mathrm{Im}\,T_{nn}(\mathbf{k},\varepsilon_{n\mathbf{k}}+i\eta)$"
             rf"   ($\eta$={ETA*RY*1e3:.0f} meV, $n_d$=1%)",y=1.00)
plt.tight_layout(rect=[0,0,1,0.89]); plt.savefig("p7_scatrate.png",dpi=140,bbox_inches="tight")
print("wrote p7_scatrate.png",flush=True)

for sd,edge,lbl in [("V",VBMe,"VBM"),("C",CBMe,"CBM")]:
    print(f"--- {lbl} side (hbar/tau in meV @ n_d=1%; tau=hbar/rate) ---")
    for d0,d1 in [(0.0,0.05),(0.05,0.15),(0.15,0.30)]:
        m=(side==sd)&(np.abs(eps-edge)>=d0)&(np.abs(eps-edge)<d1)
        if m.sum()==0: continue
        a,b,c=rM[m].mean(),rV[m].mean(),rB[m].mean()
        print(f"  |e-edge| in [{d0:.2f},{d1:.2f}) eV  ({m.sum():4d} states): "
              f"Born |M|^2 {c:8.3f} meV (tau {HBAR/c:8.1f} fs) | full-T M {a:8.3f} meV "
              f"(tau {HBAR/a:8.1f} fs) | full-T V~ {b:8.3f} meV (tau {HBAR/b:8.1f} fs) "
              f"| V~/M {b/a:6.3f} | M/Born {a/c:6.3f}")
