#!/usr/bin/env python3
"""FULL-ORDER static-omega0 rest dressing by DIRECT resolvent inversion on the explicit-60 data.

P = the 11-band window (NSCF bands 7-17, the previous Wannier/active space, N_P=1584).
Q = the remaining explicit bands (18-70, ~50/k, N_Q=7174).
  Sigma_full(w) = W_PQ (w - E_Q - W_QQ)^{-1} W_QP,   W = M/N_k (physical, eV)
via ONE eigendecomposition  E_Q + W_QQ = U diag(lam) U^H,  B = U^H W_QP:
  Sigma(w) = B^H diag(1/(w-lam)) B          -> every w costs one GEMM.
No Krylov, no stagnation; the truncated rest (to +21 eV) is band-converged to ~15 meV
(explicit-60 study).  Outputs: ladder ratios r3/r4/r5 (bare-G0 series), the QHQ spectrum
near the gap (MINRES-stagnation diagnostic), defect levels per treatment, self-consistent
e level, and vtilde_block_fesh60.dat (same format as vtilde_block.dat: M raw, Sgblk =
N_k*Sigma_phys at omega0, Vtilde = M + Sgblk) for downstream spectral-function use."""
import numpy as np
from scipy.io import FortranFile
RY=13.605693122994; OM0=-5.955; MID=-5.1

fb=FortranFile("vtilde_explicit60.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); om0r,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:2*N_A]
eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()   # raw Ry
eps=eact*RY; Nk=float(nk)
iP=np.where(a2band<=17)[0]; iQ=np.where(a2band>=18)[0]
NP_,NQ_=len(iP),len(iQ)
print(f"N_A={N_A}  P(7-17)={NP_}  Q(18-70)={NQ_}",flush=True)
VBM=eps[(a2band<=13)]; VBM=VBM[VBM<MID].max(); CBM=eps[(a2band>=14)&(eps>MID)].min()
print(f"VBM={VBM:.3f} CBM={CBM:.3f} gap={CBM-VBM:.3f} eV",flush=True)

W=(M*RY)/Nk
WPP=W[np.ix_(iP,iP)]; WPQ=W[np.ix_(iP,iQ)]; WQP=W[np.ix_(iQ,iP)]; WQQ=W[np.ix_(iQ,iQ)]
epsP=eps[iP]; epsQ=eps[iQ]

# sanity: explicit-60 P-block vs the archived 11-band block M
try:
    fb=FortranFile("vtilde_block.dat","r")
    nab,_,_,_=fb.read_ints(np.int32); fb.read_reals(np.float64); fb.read_ints(np.int32); fb.read_reals(np.float64)
    Mref=fb.read_record(np.complex128).reshape((nab,nab),order="F"); fb.close()
    print(f"sanity M_PP vs vtilde_block.dat: max|diff| = {np.abs(M[np.ix_(iP,iP)]-Mref).max():.2e} Ry",flush=True)
except Exception as ex:
    print("(sanity vs vtilde_block.dat skipped:",ex,")",flush=True)

# ---- ladder ratios with the BARE propagator (series divergence assessment, Q to band 70) ----
g0=1.0/(OM0-epsQ)
X=g0[:,None]*WQP; S2=WPQ@X
X=g0[:,None]*(WQQ@X); S3=WPQ@X
X=g0[:,None]*(WQQ@X); S4=WPQ@X
X=g0[:,None]*(WQQ@X); S5=WPQ@X
fro=lambda A: float(np.linalg.norm(A))
print(f"ladder (bare G0, omega0): ||S2..S5|| = {fro(S2):.3f} {fro(S3):.3f} {fro(S4):.3f} {fro(S5):.3f} eV")
print(f"  r3={fro(S3)/fro(S2):.3f}  r4={fro(S4)/fro(S3):.3f}  r5={fro(S5)/fro(S4):.3f}   (per-order ratio -> rho)",flush=True)

# ---- ONE eigendecomposition of the rest Hamiltonian H_QQ = E_Q + W_QQ ----
print("diagonalizing H_QQ (7174^2 zheevd, single thread ~15-25 min) ...",flush=True)
HQQ=np.diag(epsQ.astype(complex))+WQQ; HQQ=(HQQ+HQQ.conj().T)/2
lam,U=np.linalg.eigh(HQQ)
B=U.conj().T@WQP
# stagnation diagnostic: QHQ spectrum near the gap / omega0
near=lam[(lam>VBM-1.0)&(lam<CBM+1.0)]
print(f"QHQ spectrum: min(lam)={lam.min():.3f}  eigenvalues in [VBM-1,CBM+1] = {len(near)}")
if len(near): print("  nearest to omega0:", ", ".join(f"{x:.3f}" for x in near[np.argsort(np.abs(near-OM0))[:6]]))
print(f"  min|lam - omega0| = {np.abs(lam-OM0).min():.4f} eV   (A=Q(H0+dV-w0)Q near-zero mode scale)",flush=True)

def Sigma(w):
    SB=B*(1.0/(w-lam))[:,None]
    return B.conj().T@SB

def ingap(Sig):
    H=np.diag(epsP.astype(complex))+WPP+Sig; H=(H+H.conj().T)/2
    ww=np.linalg.eigvalsh(H)
    return ww[(ww>VBM+3e-3)&(ww<CBM-3e-3)]
fmt=lambda ig: ", ".join(f"{x:.3f}({x-VBM:+.3f})" for x in ig) if len(ig) else "none"

Z=np.zeros_like(WPP)
Sful=Sigma(OM0)
print(f"\nSigma_full(omega0) Hermiticity: max|S-S^H| = {np.abs(Sful-Sful.conj().T).max():.2e} eV",flush=True)
print("\n=== defect levels (eV, rel VBM), P=11 bands, Q=18-70, static omega0=-5.955 ===")
print(f"  bare M               : {fmt(ingap(Z))}")
print(f"  + Sigma2(w0)         : {fmt(ingap(S2))}")
print(f"  + Sigma2+Sigma3(w0)  : {fmt(ingap(S2+S3))}")
print(f"  + Sigma_FULL(w0)     : {fmt(ingap(Sful))}",flush=True)

# self-consistent omega for the e level
om=-4.7
for it in range(40):
    ig=ingap(Sigma(om))
    if len(ig)==0: print("  (self-consistent: no in-gap level at om=%.3f)"%om); break
    e=ig[np.argmin(np.abs(ig-om))]
    if abs(e-om)<1e-5: om=e; break
    om=0.5*om+0.5*e
print(f"  self-consistent e    : {om:.3f}({om-VBM:+.3f})")
print(f"  references: explicit-60 all-band +1.205 | explicit-21 +1.348 | DFT +1.19 | old 2nd-order(full rest) +0.36",flush=True)

# ---- write the dressed block in vtilde_block.dat format (for spectral functions etc) ----
Sg_raw=(Nk*Sful/RY).astype(np.complex128)          # file convention: Sgblk = N_k * Sigma_phys (Ry)
Mp=M[np.ix_(iP,iP)].astype(np.complex128)
Vt=Mp+Sg_raw
fo=FortranFile("vtilde_block_fesh60.dat","w")
fo.write_record(np.array([NP_,nk,nk_use,nbndskip],np.int32))
fo.write_record(np.array([OM0/RY,wmin,wmax],np.float64))
fo.write_record(np.concatenate([a2k[iP],a2band[iP]]).astype(np.int32))
fo.write_record(eact[iP].astype(np.float64))
fo.write_record(Mp.reshape(-1,order="F"))
fo.write_record(Sg_raw.reshape(-1,order="F"))
fo.write_record(Vt.reshape(-1,order="F"))
fo.close()
print("\nwrote vtilde_block_fesh60.dat  (M raw, Sgblk=N_k*Sigma_full(omega0), Vtilde=M+Sgblk; 1584^2)",flush=True)
np.savez("fesh60_partialT.npz",lam=lam,VBM=VBM,CBM=CBM,epsP=epsP)
print("done")
