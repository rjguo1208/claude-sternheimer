# ===========================================================================
#  EDT — Electron-Defect T-matrix.  Builds edt.x by reusing the compiled EDI
#  objects + QE libraries.  Source lives in the claude-sternheimer repo; the
#  QE tree (with EDI already built) is referenced by absolute path.
#  Compile with the Anvil toolchain first:
#     module reset
#     module load aocc/3.1.0 openmpi/4.1.6
#     module load amdblis/3.0 amdlibflame/3.0 amdlibm/3.0 fftw
# ===========================================================================
QEROOT  = /anvil/projects/x-che190065/rjguo/qe-7.5
EDIROOT = $(QEROOT)/edi-dev

include $(QEROOT)/make.inc

# module search path: QE base modules + PW + the compiled EDI modules
MODFLAGS = $(BASEMOD_FLAGS) $(MOD_FLAG)$(QEROOT)/PW/src $(MOD_FLAG)$(EDIROOT)/src

# QE libraries (LR_Modules added vs EDI, for the Phase-2 Sternheimer solver)
PWOBJS = $(QEROOT)/PW/src/libpw.a
LRMOD  = $(QEROOT)/LR_Modules/liblrmod.a
QEMODS = $(QEROOT)/Modules/libqemod.a \
         $(QEROOT)/upflib/libupf.a \
         $(QEROOT)/KS_Solvers/libks_solvers.a \
         $(QEROOT)/FFTXlib/src/libqefft.a \
         $(QEROOT)/UtilXlib/libutil.a \
         $(QEROOT)/XClib/xc_lib.a \
         $(QEROOT)/LAXlib/libqela.a \
         $(QEROOT)/MBD/libmbd.a \
         $(QEROOT)/dft-d3/libdftd3qe.a
W90LIB = $(QEROOT)/external/wannier90/libwannier.a

# resolve circular deps among the static QE/LR/PW archives
LIBGRP = -Wl,--start-group $(PWOBJS) $(LRMOD) $(QEMODS) $(W90LIB) -Wl,--end-group

# compiled EDI objects we reuse (everything edi.x links except edi.o itself)
EDIOBJ = $(EDIROOT)/src/edic_mod.o $(EDIROOT)/src/delta_weights.o $(EDIROOT)/src/bz_symmetry.o \
         $(EDIROOT)/src/wigner_seitz.o $(EDIROOT)/src/wan2bloch.o $(EDIROOT)/src/wannierize.o \
         $(EDIROOT)/src/edi_read_hr.o $(EDIROOT)/src/ep_constants.o $(EDIROOT)/src/io_var_edi.o \
         $(EDIROOT)/src/input_edi.o $(EDIROOT)/src/wann_common.o $(EDIROOT)/src/global_var.o \
         $(EDIROOT)/src/low_lvl_edi.o $(EDIROOT)/src/parallelism_edi.o $(EDIROOT)/src/openfilepw_edi.o \
         $(EDIROOT)/src/io_edi.o $(EDIROOT)/src/kfold_edi.o $(EDIROOT)/src/edi_input.o \
         $(EDIROOT)/src/edbloch2wan.o $(EDIROOT)/src/range_sep.o $(EDIROOT)/src/lr_matelem.o \
         $(EDIROOT)/src/ed_coarse.o $(EDIROOT)/src/transport.o $(EDIROOT)/src/tt_compress.o \
         $(EDIROOT)/src/get_betavkb.o $(EDIROOT)/src/get_vloc_onthefly.o \
         $(EDIROOT)/src/pw2wan_edi.o $(EDIROOT)/src/edi_pw2wan.o

EDT_OBJ = edt_input.o edt_partition.o edt_wannier.o edt_source.o edt.o

.DEFAULT_GOAL := edt.x

%.o : %.f90
	$(MPIF90) $(F90FLAGS) -c $< -o $@

edt.x : $(EDT_OBJ)
	$(LD) $(LDFLAGS) -o edt.x $(EDT_OBJ) $(EDIOBJ) $(LIBGRP) $(LIBOBJS) $(QELIBS)

# intra-package module dependencies
edt_partition.o : edt_partition.f90
edt_wannier.o   : edt_wannier.f90
edt_source.o    : edt_source.f90
edt.o           : edt.f90 edt_input.o edt_partition.o edt_wannier.o edt_source.o

clean:
	rm -f *.o *.mod edt.x
