#!/usr/bin/env python3
"""Diagonal of the FULL T-matrix at the VBM state vs the active-summation energy omega.

  T_PP(omega) = [1 - V~ G^A(omega)]^{-1} V~ ,  G^A_a = (1/N_k)/(omega - eps_a + i*eta)

Diagonal element T_PP(K,K;omega), VBM state a_K=(k=K, band 13) -- i.e. initial k = final k.

Closed form (exact): with z=omega+i*eta and the omega-INDEPENDENT Hermitian effective Hamiltonian
H_eff = diag(eps_a) + (1/N_k) V~  [eigenpairs (lambda_n, S)],
    T_PP(K,K;omega) = (z - eps_K) * sum_n c_n/(z - lambda_n),   c_n = S[K,n]*(S[:,n]^dag . V~[:,K]).
So one eigendecomposition gives every omega as an O(N) sum. The poles lambda_n ARE the T-matrix
resonances (the active manifold dressed by V~). Validated against a direct solve at omega0."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; ETA_EV=0.05

fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]
eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.read_record(np.complex128)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F"); fw.close()
Ki=int(np.argmin([min(((xk[0,k]-1/3+p)**2+(xk[1,k]-1/3+q)**2) for p in(-1,0,1) for q in(-1,0,1)) for k in range(nk)]))
slot=int(np.where(a2band.reshape(nk,nb)[Ki]==13)[0][0]); aK=Ki*nb+slot   # band 13 = VBM (was 17 = conduction)
eta=ETA_EV/RY
print(f"VBM state: K at k-index {Ki}, band 13, a={aK}; eps_aK={eact[aK]*RY:.4f} eV  omega0={omega0*RY:.4f} eV")
print(f"bare diagonals: |V~_KK|={abs(V[aK,aK]):.4f}  |M_KK|={abs(M[aK,aK]):.4f} Ry")

def diag_curve(Mat, oms):
    lam,S=np.linalg.eigh(np.diag(eact)+Mat/nk)            # H_eff = E + (1/Nk) Mat (Hermitian)
    c=S[aK,:]*(S.conj().T@Mat[:,aK])                       # c_n = S[K,n]*(S[:,n]^dag . Mat[:,K])
    z=oms+1j*eta
    return (z-eact[aK])*((c[None,:]/(z[:,None]-lam[None,:])).sum(1)), lam

oms=np.linspace(eact.min()-0.3/RY, omega0+1.0/RY, 2000)
Tpp,lamV=diag_curve(V,oms); Tm,lamM=diag_curve(M,oms); omv=oms*RY

# --- validate the closed form against a direct solve at omega0 ---
g=(1.0/nk)/(omega0-eact+1j*eta)
Tdir=np.linalg.solve(np.eye(N_A)-V*g[None,:], V[:,aK])[aK]
Tcf=diag_curve(V,np.array([omega0]))[0][0]
print(f"validate at omega0: direct solve T_PP(K,K)={Tdir:.5f} ; closed form={Tcf:.5f} ; |diff|={abs(Tdir-Tcf):.2e}")
iv=int(np.argmin(np.abs(oms-omega0)))
print(f"at omega0(VBM): T_PP(K,K)={Tpp[iv].real:+.4f}{Tpp[iv].imag:+.4f}i  |T_PP|={abs(Tpp[iv]):.4f}  |T_M|={abs(Tm[iv]):.4f} Ry")
near=lamV[(lamV>=(omega0-1.5/RY))&(lamV<=(omega0+0.3/RY))]
print(f"resonance poles lambda_n near VBM (eV): {np.round(near*RY,4)}")

fig,ax=plt.subplots(1,2,figsize=(13.5,5.2))
ax[0].axvspan(wmin*RY,wmax*RY,color="0.92",label="active window")
ax[0].plot(omv,Tpp.real,c="#1f77b4",lw=1.3,label=r"Re $T_{PP}$")
ax[0].plot(omv,Tpp.imag,c="#d62728",lw=1.3,label=r"Im $T_{PP}$")
ax[0].plot(omv,np.abs(Tpp),c="k",lw=2.0,label=r"$|T_{PP}|$")
ax[0].axvline(omega0*RY,ls="--",c="green",lw=1.5,label=r"$\omega_0=$VBM")
ax[0].axhline(abs(V[aK,aK]),ls=":",c="gray",label=r"bare $|\tilde V_{KK}|$ ($\omega\!\to\!\infty$)")
ax[0].set_xlabel(r"$\omega$  (eV)"); ax[0].set_ylabel(r"$T_{PP}(K,K;\omega)$  (Ry)")
ax[0].set_title("Full $T$-matrix diagonal at the VBM state vs active-summation $\\omega$")
ax[0].legend(fontsize=8,loc="upper left"); ax[0].grid(alpha=.3)
ax[1].axhline(0,c="0.6",lw=0.6,zorder=0)
ax[1].plot(omv,Tpp.real,c="#1f77b4",lw=1.7,label=r"Re $T_{PP}$")
ax[1].plot(omv,Tpp.imag,c="#d62728",lw=1.7,label=r"Im $T_{PP}$")
ax[1].plot(omv,np.abs(Tpp),c="k",lw=2.1,label=r"$|T_{PP}|$")
ax[1].plot(omv,Tm.real,c="#1f77b4",lw=1.4,ls="--",label=r"Re $T_M$")
ax[1].plot(omv,Tm.imag,c="#d62728",lw=1.4,ls="--",label=r"Im $T_M$")
ax[1].plot(omv,np.abs(Tm),c="k",lw=1.4,ls="--",label=r"$|T_M|$")
ax[1].axvline(omega0*RY,ls="--",c="green",lw=1.3,label=r"$\omega_0=$VBM")
ax[1].set_xlim(omega0*RY-1.5, omega0*RY+0.6)
ax[1].set_xlabel(r"$\omega$  (eV)"); ax[1].set_ylabel(r"$T(K,K;\omega)$  (Ry)")
ax[1].set_title("Band-edge zoom — solid: with rest-space $T_{PP}$, dashed: no rest-space $T_M$")
ax[1].legend(fontsize=7,ncol=2); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig("p6_tdiag_omega.png",dpi=130); print("wrote p6_tdiag_omega.png")
