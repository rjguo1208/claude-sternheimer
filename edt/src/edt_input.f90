MODULE edt_input
  !---------------------------------------------------------------------------
  !  EDT — Electron-Defect T-matrix.  Namelist &edt_nml.
  !  Extends EDI's input with the active/rest + Sternheimer T-matrix controls
  !  (plan.md §12).  All keys are read on ionode and broadcast to all ranks.
  !---------------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  SAVE
  PUBLIC

  ! ---- primitive cell (QE NSCF save) ----
  CHARACTER(LEN=256) :: edi_prefix = 'pwscf'
  CHARACTER(LEN=256) :: edi_outdir = './'

  ! ---- coarse k-grid (must match the primitive NSCF) ----
  INTEGER :: coarse_nk1 = 0, coarse_nk2 = 0, coarse_nk3 = 0

  ! ---- supercell difference potential (cube files from extract_pot.x) ----
  CHARACTER(LEN=256) :: potfile_d = ' '
  CHARACTER(LEN=256) :: potfile_p = ' '
  CHARACTER(LEN=256) :: pot_align = 'vacuum'        ! 'vacuum' | 'core' | 'none'
  REAL(dp) :: defect_center(3) = 0.0_dp
  REAL(dp) :: core_align_radius = 2.0_dp

  ! ---- Wannier data ----
  INTEGER :: nbndsub = 0                            ! number of Wannier functions (active manifold)
  CHARACTER(LEN=256) :: filukk   = 'filukk'         ! combined rotation U(k)
  CHARACTER(LEN=256) :: hr_seedname = ' '           ! seedname for <seed>_hr.dat (default = edi_prefix)

  ! ---- range separation (long-range Coulomb tail; route SR through Sternheimer) ----
  LOGICAL  :: range_sep   = .FALSE.
  CHARACTER(LEN=256) :: rhofile_d = ' '
  CHARACTER(LEN=256) :: rhofile_p = ' '
  LOGICAL  :: coulomb_2d  = .TRUE.
  REAL(dp) :: alpha_gauss = 0.0_dp

  ! ---- T-matrix / Sternheimer controls (NEW; plan.md §12) ----
  LOGICAL  :: do_tmatrix      = .TRUE.
  REAL(dp) :: active_win_min  = -9999.0_dp          ! active window [eV] defining A(k)
  REAL(dp) :: active_win_max  =  9999.0_dp
  REAL(dp) :: omega0          =  9999.0_dp          ! sentinel: default = VBM of the active manifold
  REAL(dp) :: eta             =  0.01_dp            ! complex Sternheimer shift omega0 + i*eta [eV]
  INTEGER  :: rest_nk1 = 0, rest_nk2 = 0, rest_nk3 = 0   ! full-BZ rest grid (0 -> use coarse grid)
  REAL(dp) :: sternheimer_thr = 1.0d-10
  INTEGER  :: dress_order     = 0                   ! V_QQ ladder rungs (0 = coupling-2nd-order)
  REAL(dp) :: dress_tol       = 1.0d-3
  LOGICAL  :: active_resum    = .TRUE.
  CHARACTER(LEN=16) :: resum_grid = 'coarse'        ! 'coarse' (B1) | 'fine' (B2 Dyson)
  CHARACTER(LEN=16) :: rest_split = 'complex'       ! 'complex' (ccgsolve) | 'pm' (split R^+/-)

  ! ---- P3 full off-diagonal block assembly (MPI, pool-parallel; plan.md §P3) ----
  LOGICAL  :: do_full_block    = .FALSE.            ! assemble the full Vtilde block on the active manifold
  INTEGER  :: block_nk         = 0                  ! restrict manifold+channels to k=1..block_nk (0 = full BZ)
  INTEGER  :: block_single_band= 0                  ! validation: single-ket diag check at (band, k); 0 = full block
  INTEGER  :: block_single_ki  = 1
  CHARACTER(LEN=256) :: vtilde_outfile = 'vtilde_block.dat'
  LOGICAL  :: dump_wann        = .FALSE.            ! P5-b: dump xk_cryst + U(k) to wann_data.dat (run in audit mode)

  NAMELIST / edt_nml / &
       edi_prefix, edi_outdir, coarse_nk1, coarse_nk2, coarse_nk3, &
       potfile_d, potfile_p, pot_align, defect_center, core_align_radius, &
       nbndsub, filukk, hr_seedname, &
       range_sep, rhofile_d, rhofile_p, coulomb_2d, alpha_gauss, &
       do_tmatrix, active_win_min, active_win_max, omega0, eta, &
       rest_nk1, rest_nk2, rest_nk3, sternheimer_thr, dress_order, dress_tol, &
       active_resum, resum_grid, rest_split, &
       do_full_block, block_nk, block_single_band, block_single_ki, vtilde_outfile, dump_wann

CONTAINS

  SUBROUTINE read_edt_input()
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : world_comm
    IMPLICIT NONE
    INTEGER :: ios

    ios = 0
    IF (ionode) THEN
       CALL input_from_file()
       READ(5, NML=edt_nml, IOSTAT=ios)
       IF (ios < 0) ios = 0      ! EOF after the namelist is fine
    ENDIF
    CALL mp_bcast(ios, ionode_id, world_comm)
    IF (ios > 0) CALL errore('read_edt_input', 'error reading &edt_nml namelist', ABS(ios))

    ! default the hr seedname to the prefix
    IF (ionode .AND. LEN_TRIM(hr_seedname) == 0) hr_seedname = edi_prefix

    CALL mp_bcast(edi_prefix,       ionode_id, world_comm)
    CALL mp_bcast(edi_outdir,       ionode_id, world_comm)
    CALL mp_bcast(coarse_nk1,       ionode_id, world_comm)
    CALL mp_bcast(coarse_nk2,       ionode_id, world_comm)
    CALL mp_bcast(coarse_nk3,       ionode_id, world_comm)
    CALL mp_bcast(potfile_d,        ionode_id, world_comm)
    CALL mp_bcast(potfile_p,        ionode_id, world_comm)
    CALL mp_bcast(pot_align,        ionode_id, world_comm)
    CALL mp_bcast(defect_center,    ionode_id, world_comm)
    CALL mp_bcast(core_align_radius,ionode_id, world_comm)
    CALL mp_bcast(nbndsub,          ionode_id, world_comm)
    CALL mp_bcast(filukk,           ionode_id, world_comm)
    CALL mp_bcast(hr_seedname,      ionode_id, world_comm)
    CALL mp_bcast(range_sep,        ionode_id, world_comm)
    CALL mp_bcast(rhofile_d,        ionode_id, world_comm)
    CALL mp_bcast(rhofile_p,        ionode_id, world_comm)
    CALL mp_bcast(coulomb_2d,       ionode_id, world_comm)
    CALL mp_bcast(alpha_gauss,      ionode_id, world_comm)
    CALL mp_bcast(do_tmatrix,       ionode_id, world_comm)
    CALL mp_bcast(active_win_min,   ionode_id, world_comm)
    CALL mp_bcast(active_win_max,   ionode_id, world_comm)
    CALL mp_bcast(omega0,           ionode_id, world_comm)
    CALL mp_bcast(eta,              ionode_id, world_comm)
    CALL mp_bcast(rest_nk1,         ionode_id, world_comm)
    CALL mp_bcast(rest_nk2,         ionode_id, world_comm)
    CALL mp_bcast(rest_nk3,         ionode_id, world_comm)
    CALL mp_bcast(sternheimer_thr,  ionode_id, world_comm)
    CALL mp_bcast(dress_order,      ionode_id, world_comm)
    CALL mp_bcast(dress_tol,        ionode_id, world_comm)
    CALL mp_bcast(active_resum,     ionode_id, world_comm)
    CALL mp_bcast(resum_grid,       ionode_id, world_comm)
    CALL mp_bcast(rest_split,       ionode_id, world_comm)
    CALL mp_bcast(do_full_block,    ionode_id, world_comm)
    CALL mp_bcast(block_nk,         ionode_id, world_comm)
    CALL mp_bcast(block_single_band,ionode_id, world_comm)
    CALL mp_bcast(block_single_ki,  ionode_id, world_comm)
    CALL mp_bcast(vtilde_outfile,   ionode_id, world_comm)
    CALL mp_bcast(dump_wann,        ionode_id, world_comm)

    IF (ionode) THEN
       WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
       WRITE(stdout,'(5X,A)')   'EDT input (&edt_nml)'
       WRITE(stdout,'(5X,A)')   REPEAT('=',64)
       WRITE(stdout,'(5X,A,A)')      'edi_prefix       = ', TRIM(edi_prefix)
       WRITE(stdout,'(5X,A,A)')      'edi_outdir       = ', TRIM(edi_outdir)
       WRITE(stdout,'(5X,A,3I5)')    'coarse_nk        = ', coarse_nk1, coarse_nk2, coarse_nk3
       WRITE(stdout,'(5X,A,I5)')     'nbndsub (N_W)    = ', nbndsub
       WRITE(stdout,'(5X,A,A)')      'filukk           = ', TRIM(filukk)
       WRITE(stdout,'(5X,A,A)')      'hr_seedname      = ', TRIM(hr_seedname)
       WRITE(stdout,'(5X,A,2F10.4)') 'active window eV = ', active_win_min, active_win_max
       IF (omega0 > 9000.0_dp) THEN
          WRITE(stdout,'(5X,A)')     'omega0           = (default: VBM of active manifold)'
       ELSE
          WRITE(stdout,'(5X,A,F10.4)') 'omega0 eV        = ', omega0
       ENDIF
       WRITE(stdout,'(5X,A,F10.5)')  'eta eV           = ', eta
       WRITE(stdout,'(5X,A,3I5)')    'rest_nk          = ', rest_nk1, rest_nk2, rest_nk3
       WRITE(stdout,'(5X,A,L2)')     'do_tmatrix       = ', do_tmatrix
       WRITE(stdout,'(5X,A,I3,A,L2)')'dress_order      = ', dress_order, '   active_resum = ', active_resum
       WRITE(stdout,'(5X,A,A,A,A)')  'resum_grid       = ', TRIM(resum_grid), '   rest_split = ', TRIM(rest_split)
       WRITE(stdout,'(5X,A,L2,A,A)') 'range_sep        = ', range_sep, '   potfiles: ', &
            TRIM(potfile_d)//' , '//TRIM(potfile_p)
       WRITE(stdout,'(5X,A,L2,A,I5)')'do_full_block    = ', do_full_block, '   block_nk = ', block_nk
       IF (block_single_band > 0) &
          WRITE(stdout,'(5X,A,I5,A,I5)') 'block single-ket band = ', block_single_band, '  ki = ', block_single_ki
       WRITE(stdout,'(5X,A)')   REPEAT('=',64)
       FLUSH(stdout)
    ENDIF
  END SUBROUTINE read_edt_input

END MODULE edt_input
