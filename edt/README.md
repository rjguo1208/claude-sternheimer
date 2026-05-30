# EDT — Electron-Defect T-matrix (Phase 0)

QE plug-in implementing the downfolding + Sternheimer electron–defect *T*-matrix
(see `../plan.md`). Phase 0 = scaffold, build, Stage-A inputs, and the active/rest
partition audit. **No Sternheimer solve yet** (Phase 2+).

## Layout
```
edt/src/   edt_input.f90      &edt_nml namelist (read + broadcast)
           edt_wannier.f90    minimal robust filukk reader (u_kc rotation)
           edt_partition.f90  active/rest split, apply_Qproj, audit print
           edt.f90            main: read NSCF + filukk + H_W(R), partition, audit, load ΔV
           makefile           links against compiled EDI objects + QE libs (absolute QEROOT)
edt/run/   edt_audit.in       audit only (fast)
           edt_full.in        audit + supercell ΔV cube load
```

It **reuses the compiled EDI objects** in `$(EDIROOT)/src/*.o` and the QE libraries
(adds `LR_Modules/liblrmod.a` for the Phase-2 Sternheimer solver). Source lives here
in the repo; `QEROOT`/`EDIROOT` in the makefile point at the QE tree.

## Build (Anvil)
```bash
module reset
module load aocc/3.1.0 openmpi/4.1.6
module load amdblis/3.0 amdlibflame/3.0 amdlibm/3.0 fftw
cd edt/src && make
```
(Load modules with plain redirection, not a pipe — piping `module load` into `tail`
runs it in a subshell and drops `LIBRARY_PATH`, breaking the `-lfftw3` link.)

## Run (login node, single rank is fine for the audit)
```bash
cd edt/run
mpirun -np 1 ../src/edt.x -i edt_audit.in > edt_audit.out
```

## Data (MoS₂, high-band NSCF for the rest space)
- primitive NSCF (150 bands, 12×12): `…/T-matrix/data/primitive_highbands_12x12/dout/mos2.save`
- Wannier rotation / Hamiltonian: `…/T-matrix/data/filukk`, `…/T-matrix/data/mos2_hr.dat`
- supercell ΔV cubes: `…/T-matrix/V_d.cube`, `…/T-matrix/V_p.cube`

The Wannierization filukk was written with nbnd=17 (6 semicore + 11 valence); the NSCF
here has nbnd=150. EDI's `read_filukk_edi` assumes the trailing `excluded_band` block has
the *current* nbnd and overruns, so `edt_wannier.f90` reads only `nbndep`/`ibndkept`/`u_kc`.

Partition for this system: **active = 11 Wannier valence bands (7–17); rest = 133 bands
(18–150)**; 6 semicore (1–6) dropped. Wannier interpolation matches NSCF to ~2×10⁻⁵ eV.
