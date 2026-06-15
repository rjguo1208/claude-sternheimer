#!/usr/bin/env python3
"""Parametrized full-order static-omega0 rest dressing (direct resolvent) for ANY
explicit-60 raw-M file -> dressed 11-band block for spectral functions.
Same math as fesh60_partialT.py (P=NSCF bands 7-17, Q=18-70, Sigma=W_PQ(w-E_Q-W_QQ)^-1 W_QP
via one H_QQ eigendecomposition); only the input/output filenames are arguments.

usage:  fesh60_partialT_sub.py <in_explicit60.dat> <out_block.dat> [TAG] [OM0_eV]
"""
import sys, numpy as np
from scipy.io import FortranFile
RY=13.605693122994; MID=-5.1
FN  = sys.argv[1] if len(sys.argv)>1 else "vtilde_explicit60.dat"
OUT = sys.argv[2] if len(sys.argv)>2 else "vtilde_block_fesh60.dat"
TAG = sys.argv[3] if len(sys.argv)>3 else "sub"
OM0 = float(sys.argv[4]) if len(sys.argv)>4 else -5.955

fb=FortranFile(FN,"r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); om0r,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:2*N_A]
eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
eps=eact*RY; Nk=float(nk)
iP=np.where(a2band<=17)[0]; iQ=np.where(a2band>=18)[0]; NP_,NQ_=len(iP),len(iQ)
print(f"[{TAG}] N_A={N_A}  P(7-17)={NP_}  Q(18-70)={NQ_}  OM0={OM0}",flush=True)
VBM=eps[(a2band<=13)]; VBM=VBM[VBM<MID].max(); CBM=eps[(a2band>=14)&(eps>MID)].min()
print(f"[{TAG}] VBM={VBM:.3f} CBM={CBM:.3f} gap={CBM-VBM:.3f} eV",flush=True)
W=(M*RY)/Nk
WPP=W[np.ix_(iP,iP)]; WPQ=W[np.ix_(iP,iQ)]; WQP=W[np.ix_(iQ,iP)]; WQQ=W[np.ix_(iQ,iQ)]
epsP=eps[iP]; epsQ=eps[iQ]

print(f"[{TAG}] diagonalizing H_QQ ({NQ_}^2 zheevd) ...",flush=True)
HQQ=np.diag(epsQ.astype(complex))+WQQ; HQQ=(HQQ+HQQ.conj().T)/2
lam,U=np.linalg.eigh(HQQ); B=U.conj().T@WQP
def Sigma(w): return B.conj().T@(B*(1.0/(w-lam))[:,None])
def ingap(Sig):
    H=np.diag(epsP.astype(complex))+WPP+Sig; H=(H+H.conj().T)/2
    ww=np.linalg.eigvalsh(H); return ww[(ww>VBM+3e-3)&(ww<CBM-3e-3)]
fmt=lambda ig: ", ".join(f"{x:.3f}({x-VBM:+.3f})" for x in ig) if len(ig) else "none"
Z=np.zeros_like(WPP); Sful=Sigma(OM0)
print(f"[{TAG}] Sigma_full(omega0) Hermiticity max|S-S^H| = {np.abs(Sful-Sful.conj().T).max():.2e} eV",flush=True)
print(f"[{TAG}] === in-gap levels (eV, rel VBM), P=11 bands ===")
print(f"[{TAG}]   bare M            : {fmt(ingap(Z))}")
print(f"[{TAG}]   + Sigma_FULL(w0)  : {fmt(ingap(Sful))}",flush=True)
# self-consistent omega (scan from mid-gap)
om=0.5*(VBM+CBM)
for it in range(40):
    ig=ingap(Sigma(om))
    if len(ig)==0: print(f"[{TAG}]   self-consistent: no in-gap level near {om:.3f}"); break
    e=ig[np.argmin(np.abs(ig-om))]
    if abs(e-om)<1e-5: om=e; break
    om=0.5*om+0.5*e
else:
    pass
if 'e' in dir() and len(ingap(Sigma(om))): print(f"[{TAG}]   self-consistent lvl: {om:.3f}({om-VBM:+.3f})",flush=True)

Sg_raw=(Nk*Sful/RY).astype(np.complex128); Mp=M[np.ix_(iP,iP)].astype(np.complex128); Vt=Mp+Sg_raw
fo=FortranFile(OUT,"w")
fo.write_record(np.array([NP_,nk,nk_use,nbndskip],np.int32))
fo.write_record(np.array([OM0/RY,wmin,wmax],np.float64))
fo.write_record(np.concatenate([a2k[iP],a2band[iP]]).astype(np.int32))
fo.write_record(eact[iP].astype(np.float64))
fo.write_record(Mp.reshape(-1,order="F"))
fo.write_record(Sg_raw.reshape(-1,order="F"))
fo.write_record(Vt.reshape(-1,order="F"))
fo.close()
print(f"[{TAG}] wrote {OUT}  (M raw, Sgblk=N_k*Sigma_full(omega0), Vtilde=M+Sgblk; {NP_}^2)",flush=True)
