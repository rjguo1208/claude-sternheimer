#!/usr/bin/env python3
"""Read edt vtilde_block.dat (Fortran unformatted) and report physics diagnostics.
Record layout (see edt_sternheimer.f90::vtilde_block_mpi):
  1: int32   N_A, nkstot, nk_use, nbndskip
  2: float64 omega0_ry, win_min_ry, win_max_ry
  3: int32   a2k(1:N_A), a2band(1:N_A)
  4: float64 eact(1:N_A)            (active-state energies, Ry)
  5: complex128 Mblk  (N_A x N_A, Fortran order)   Born
  6: complex128 Sgblk (N_A x N_A, Fortran order)   rest self-energy
  7: complex128 Vblk  (N_A x N_A, Fortran order)   Vtilde = M + Sigma (Hermitized)
"""
import sys, numpy as np
from scipy.io import FortranFile

RY = 13.605693122994

f = FortranFile(sys.argv[1] if len(sys.argv) > 1 else "vtilde_block.dat", "r")
N_A, nkstot, nk_use, nbndskip = f.read_ints(np.int32)
omega0, wmin, wmax = f.read_reals(np.float64)
idx = f.read_ints(np.int32); a2k = idx[:N_A]; a2band = idx[N_A:]
eact = f.read_reals(np.float64)
M  = f.read_record(np.complex128).reshape((N_A, N_A), order="F")
Sg = f.read_record(np.complex128).reshape((N_A, N_A), order="F")
V  = f.read_record(np.complex128).reshape((N_A, N_A), order="F")
f.close()

dM, dS, dV = np.diag(M).real, np.diag(Sg).real, np.diag(V).real
print(f"N_A={N_A}  nkstot={nkstot}  nk_use={nk_use}  nbndskip={nbndskip}")
print(f"omega0={omega0*RY:.4f} eV   window=[{wmin*RY:.3f},{wmax*RY:.3f}] eV")
print(f"active states: {N_A}  ({N_A//nk_use} bands x {nk_use} k)")

print("\n--- Hermiticity (from file) ---")
print(f"  max|V - V^H| = {np.abs(V - V.conj().T).max():.2e} Ry   max|V| = {np.abs(V).max():.3e} Ry")

print("\n--- magnitudes (Frobenius norms, Ry) ---")
nM, nS, nV = np.linalg.norm(M), np.linalg.norm(Sg), np.linalg.norm(V)
print(f"  ||M||_F={nM:.3f}   ||Sigma||_F={nS:.3f}   ||Vtilde||_F={nV:.3f}")
print(f"  rest-dressing strength ||Sigma||/||M|| = {nS/nM:.3f}")
print(f"  Vtilde vs Born change  ||V-M||/||M||   = {np.linalg.norm(V-M)/nM:.3f}")

print("\n--- diagonal V_nn (Ry) ---")
print(f"  range [{dV.min():.4f}, {dV.max():.4f}]   mean {dV.mean():.4f}   rms {np.sqrt((dV**2).mean()):.4f}")
nbB = int(np.sum(np.abs(dS) > np.abs(dM)))
nflip = int(np.sum(np.sign(dV) != np.sign(dM)))
print(f"  |Sigma_nn|>|M_nn| (strong beyond-Born): {nbB}/{N_A} states")
print(f"  sign(V_nn) != sign(M_nn): {nflip}/{N_A} states  (rest dressing flips the sign)")

print("\n--- off-diagonal (k-mixing) strength ---")
offV = nV**2 - np.sum(dV**2)
diagV = np.sum(dV**2)
print(f"  ||V_offdiag||_F/||V_diag||_F = {np.sqrt(offV/diagV):.3f}   (intra+inter-k coupling vs on-site)")

print("\n--- Hermitian eigenvalue spectra (Ry) — physical scattering strengths ---")
wM = np.linalg.eigvalsh((M+M.conj().T)/2)
wV = np.linalg.eigvalsh(V)
print(f"  Born   M: [{wM.min():.3f}, {wM.max():.3f}]   (5 largest |.|: {np.sort(np.abs(wM))[-5:][::-1].round(3)})")
print(f"  Vtilde  : [{wV.min():.3f}, {wV.max():.3f}]   (5 largest |.|: {np.sort(np.abs(wV))[-5:][::-1].round(3)})")
print(f"  spectral radius: Born {np.abs(wM).max():.3f} Ry -> Vtilde {np.abs(wV).max():.3f} Ry")
