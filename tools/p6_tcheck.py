#!/usr/bin/env python3
"""Why do T_M and T_PP look similar at the VBM diagonal although M and V~ differ a lot?
Check the diagonal T_X(K,K) AND the full-matrix norm ||T_X|| at on- vs off-resonance omega,
for X=M (Born) and X=V~ (beyond-Born), via the closed form T_X=(z-E)(z-H_eff)^-1 X."""
import numpy as np
from scipy.io import FortranFile
RY=13.605693122994; eta=0.05/RY
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]; eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.read_record(np.complex128)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F"); fw.close()
Ki=int(np.argmin([min(((xk[0,k]-1/3+p)**2+(xk[1,k]-1/3+q)**2) for p in(-1,0,1) for q in(-1,0,1)) for k in range(nk)]))
aK=Ki*nb+int(np.where(a2band.reshape(nk,nb)[Ki]==17)[0][0])

def prep(X):
    lam,S=np.linalg.eigh(np.diag(eact)+X/nk); return lam,S,S.conj().T@X   # lam, S, S†X
def Tdiag(o,P,X):  lam,S,SdX=P; z=o+1j*eta; return (z-eact[aK])*np.sum(S[aK,:]*SdX[:,aK]/(z-lam))
def Tnorm(o,P):    lam,S,SdX=P; z=o+1j*eta; return np.linalg.norm((z-eact)[:,None]*(S@(( 1/(z-lam))[:,None]*SdX)))
PM,PV=prep(M),prep(V)

print(f"VBM state a_K={aK}; omega0(VBM)={omega0*RY:.4f} eV")
print(f"bare diagonal:  |M_KK|={abs(M[aK,aK]):.4f}  |V~_KK|={abs(V[aK,aK]):.4f}  ratio V~/M={abs(V[aK,aK])/abs(M[aK,aK]):.3f}")
print(f"bare norm:      ||M|| ={np.linalg.norm(M):.1f}   ||V~||={np.linalg.norm(V):.1f}   ratio V~/M={np.linalg.norm(V)/np.linalg.norm(M):.3f}\n")
print(f"{'omega-VBM[eV]':>12} | {'|T_M(K,K)|':>11} {'|T_PP(K,K)|':>11} {'ratio':>6} | {'||T_M||':>9} {'||T_PP||':>9} {'ratio':>6}  regime")
for dE,tag in [(+1.5,"off-res (above window)"),(+0.5,"off-res (above VBM)"),(0.0,"on-shell VBM"),(-0.5,"in-window"),(-2.0,"in-window")]:
    o=omega0+dE/RY
    tdM,tdV=abs(Tdiag(o,PM,M)),abs(Tdiag(o,PV,V)); tnM,tnV=Tnorm(o,PM),Tnorm(o,PV)
    print(f"{dE:>+12.2f} | {tdM:>11.4f} {tdV:>11.4f} {tdV/tdM:>6.2f} | {tnM:>9.1f} {tnV:>9.1f} {tnV/tnM:>6.2f}  {tag}")
