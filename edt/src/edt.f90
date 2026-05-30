PROGRAM edt
  !===========================================================================
  !  EDT — Electron-Defect T-matrix (downfolding + Sternheimer).  Phase 0:
  !  scaffold, build, Stage-A inputs, and the active/rest partition audit.
  !  No Sternheimer solve yet (that is Phase 2+).  See plan.md §P0.
  !===========================================================================
  USE kinds,            ONLY : dp
  USE io_files,         ONLY : prefix, tmp_dir
  USE wvfct,            ONLY : nbnd, et
  USE klist,            ONLY : nkstot, nks, xk
  USE cell_base,        ONLY : at
  USE constants,        ONLY : rytoev
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE mp_global,        ONLY : mp_startup, mp_global_end
  USE mp_bands,         ONLY : nbgrp
  USE noncollin_module, ONLY : noncolin, npol, lspinorb
  USE environment,      ONLY : environment_start, environment_end

  USE edt_input,        ONLY : read_edt_input, edi_prefix, edi_outdir, &
                               coarse_nk1, coarse_nk2, coarse_nk3, nbndsub, filukk, &
                               hr_seedname, active_win_min, active_win_max, omega0, eta, &
                               potfile_d, potfile_p, pot_align, defect_center, core_align_radius, &
                               range_sep, rhofile_d, rhofile_p, coulomb_2d, alpha_gauss
  USE edt_partition,    ONLY : build_active_set, print_partition_audit
  USE ed_coarse,        ONLY : load_pot_from_file, &
                               build_vcolin_aligned, build_vcolin_corealign, write_vcolin_cube
  USE edt_wannier,      ONLY : edt_read_filukk
  USE edi_read_hr,      ONLY : read_hr_file
  USE edi_pw2wan,       ONLY : edi_interp_bands
  USE range_sep,        ONLY : compute_range_separation
  USE edic_mod,         ONLY : V_d, V_p, V_colin
  USE wann_common,      ONLY : n_wannier
  USE global_var,       ONLY : nbndep, nbndskip, ibndkept

  IMPLICIT NONE
  CHARACTER(LEN=256), EXTERNAL :: trimcheck

  INTEGER :: ik, ib, j, ierr, nbndsub_loc, nrr
  INTEGER, ALLOCATABLE :: irvec(:,:), ndegen(:)
  COMPLEX(dp), ALLOCATABLE :: chw(:,:,:)
  REAL(dp), ALLOCATABLE :: xk_all(:,:), xk_cryst(:,:), eig_interp(:,:)
  REAL(dp), ALLOCATABLE :: et_ev(:,:)
  LOGICAL, ALLOCATABLE :: kept(:), is_active(:,:)
  INTEGER, ALLOCATABLE :: n_act(:), n_rest_k(:)
  INTEGER :: n_active_tot, n_rest_tot, nexcl
  REAL(dp) :: vbm, cbm, gap, max_err, err, vac_shift
  REAL(dp) :: win_min, win_max, om0

  CALL mp_startup(start_images=.TRUE.)
  CALL environment_start('EDT')
  IF (nbgrp > 1) CALL errore('edt', 'band groups not supported', nbgrp)

  ! ---- input ----
  CALL read_edt_input()
  prefix  = TRIM(edi_prefix)
  tmp_dir = trimcheck(TRIM(edi_outdir))
  CALL mp_bcast(prefix,  ionode_id, world_comm)
  CALL mp_bcast(tmp_dir, ionode_id, world_comm)

  ! ---- Stage A: read primitive NSCF ----
  IF (ionode) WRITE(stdout,'(/,5X,A)') 'Reading primitive NSCF save ...'
  CALL read_file()
  IF (ionode) THEN
     WRITE(stdout,'(5X,A,I6)') 'nkstot = ', nkstot
     WRITE(stdout,'(5X,A,I6)') 'nks    = ', nks
     WRITE(stdout,'(5X,A,I6)') 'nbnd   = ', nbnd
     WRITE(stdout,'(5X,A,I6)') 'npol   = ', npol
     IF (lspinorb) THEN
        WRITE(stdout,'(5X,A)') 'spin: SOC (noncolin=T, lspinorb=T)'
     ELSE IF (noncolin) THEN
        WRITE(stdout,'(5X,A)') 'spin: noncollinear'
     ELSE
        WRITE(stdout,'(5X,A)') 'spin: collinear (scalar)'
     ENDIF
  ENDIF
  IF (coarse_nk1*coarse_nk2*coarse_nk3 /= nkstot) &
       CALL errore('edt', 'coarse_nk1*nk2*nk3 /= nkstot', nkstot)

  ! ---- crystal-coordinate k-points (for Wannier interpolation) ----
  ALLOCATE(xk_all(3, nkstot), xk_cryst(3, nkstot))
  CALL poolcollect(3, nks, xk, nkstot, xk_all)
  xk_cryst(:,:) = xk_all(:,:)
  CALL cryst_to_cart(nkstot, xk_cryst, at, -1)
  CALL mp_bcast(xk_cryst, ionode_id, world_comm)

  ! ---- Wannier rotation U(k) + H_W(R) (reuse EDI) ----
  IF (ionode) WRITE(stdout,'(/,5X,A)') 'Reading Wannier rotation (filukk) ...'
  CALL edt_read_filukk(TRIM(filukk), nkstot, nks, nbndsub)
  IF (ionode) THEN
     WRITE(stdout,'(5X,A,I6)') 'nbndep (kept bands)  = ', nbndep
     WRITE(stdout,'(5X,A,I6)') 'n_wannier            = ', n_wannier
  ENDIF
  IF (ionode) WRITE(stdout,'(5X,A)') 'Reading H_W(R) from <seed>_hr.dat ...'
  CALL read_hr_file(TRIM(hr_seedname), nbndsub_loc, nrr, ndegen, irvec, chw)
  IF (ionode) WRITE(stdout,'(5X,A,I6,A,I8)') 'nbndsub(hr) = ', nbndsub_loc, '   nrr = ', nrr

  ! ---- sanity: Wannier-interpolated bands vs NSCF (reuse EDI) ----
  ALLOCATE(eig_interp(nbndsub_loc, nkstot))
  CALL edi_interp_bands(nbndsub_loc, nrr, irvec, ndegen, chw, nkstot, xk_cryst, eig_interp)
  IF (ionode) THEN
     max_err = 0.0_dp
     DO ik = 1, nkstot
        DO ib = 1, MIN(nbndsub_loc, nbndep)
           err = ABS(eig_interp(ib, ik) - et(ibndkept(ib), ik)*rytoev)
           max_err = MAX(max_err, err)
        ENDDO
     ENDDO
     WRITE(stdout,'(5X,A,ES12.4,A)') 'Wannier interp max |E_w90 - E_nscf| = ', max_err, ' eV'
     IF (max_err < 1.0d-3) THEN
        WRITE(stdout,'(5X,A)') 'PASS: Wannier interpolation consistent.'
     ELSE
        WRITE(stdout,'(5X,A)') 'WARN: interpolation error large (check window/froz).'
     ENDIF
  ENDIF

  ! ---- build kept(nbnd) and et in eV ----
  ! kept = full band space minus the semicore bands skipped at Wannierization.
  ! The rest space R = kept \ active then spans the high/conduction bands
  ! (the reason we use a high-band NSCF, nbnd=150).
  ALLOCATE(kept(nbnd), et_ev(nbnd, nkstot))
  kept = .TRUE.
  DO ib = 1, nbndskip
     kept(ib) = .FALSE.
  ENDDO
  nexcl = nbndskip
  et_ev(:,:) = et(:, 1:nkstot) * rytoev

  ! ---- default active window from the Wannier manifold if user left sentinel ----
  win_min = active_win_min; win_max = active_win_max
  IF (win_min <= -9000.0_dp .OR. win_max >= 9000.0_dp) THEN
     win_min = MINVAL(eig_interp) - 0.05_dp
     win_max = MAXVAL(eig_interp) + 0.05_dp
     IF (ionode) WRITE(stdout,'(/,5X,A,2F10.4,A)') &
          'active window defaulted to Wannier-manifold span = [', win_min, win_max, ' ] eV'
  ENDIF

  ! ---- active/rest partition ----
  ALLOCATE(is_active(nbnd, nkstot), n_act(nkstot), n_rest_k(nkstot))
  om0 = omega0
  CALL build_active_set(nbnd, nkstot, et_ev, kept, win_min, win_max, &
                         is_active, n_act, n_rest_k, n_active_tot, n_rest_tot, &
                         vbm, cbm, om0, gap)

  IF (ionode) CALL print_partition_audit(nbnd, nkstot, nbnd-nbndskip, n_wannier, nexcl, &
                         win_min, win_max, om0, eta, n_act, n_rest_k, &
                         n_active_tot, n_rest_tot, vbm, cbm, gap)

  ! ---- Stage A: supercell difference potential (cube files) ----
  IF (LEN_TRIM(potfile_d) > 0 .AND. LEN_TRIM(potfile_p) > 0) THEN
     IF (ionode) WRITE(stdout,'(/,5X,A)') 'Loading supercell potentials (cube) ...'
     CALL load_pot_from_file(TRIM(potfile_d), V_d)
     CALL load_pot_from_file(TRIM(potfile_p), V_p)
     ALLOCATE(V_colin(V_d%nr1 * V_d%nr2 * V_d%nr3))
     SELECT CASE (TRIM(pot_align))
     CASE ('vacuum')
        CALL build_vcolin_aligned(V_d, V_p, V_colin, SIZE(V_colin), vac_shift)
     CASE ('core')
        CALL build_vcolin_corealign(V_d, V_p, V_colin, SIZE(V_colin), vac_shift, &
                                     defect_center, core_align_radius)
     CASE ('none')
        V_colin(:) = V_d%pot(:) - V_p%pot(:); vac_shift = 0.0_dp
     CASE DEFAULT
        CALL errore('edt', 'unknown pot_align: '//TRIM(pot_align), 1)
     END SELECT
     IF (ionode) THEN
        WRITE(stdout,'(5X,A,3I6)')    'SC FFT grid          = ', V_d%nr1, V_d%nr2, V_d%nr3
        WRITE(stdout,'(5X,A,ES13.5)') 'alignment shift (Ry) = ', vac_shift
        WRITE(stdout,'(5X,A,2ES13.5)')'V_colin min/max (Ry) = ', MINVAL(V_colin), MAXVAL(V_colin)
     ENDIF

     ! optional range separation (SR/LR split)
     IF (range_sep .AND. LEN_TRIM(rhofile_d) > 0 .AND. LEN_TRIM(rhofile_p) > 0) THEN
        BLOCK
          USE edic_mod, ONLY : V_file
          TYPE(V_file) :: rho_d, rho_p
          REAL(dp), ALLOCATABLE :: V_lr(:)
          IF (ionode) WRITE(stdout,'(5X,A)') 'Range separation: loading charge-density cubes ...'
          CALL load_pot_from_file(TRIM(rhofile_d), rho_d)
          CALL load_pot_from_file(TRIM(rhofile_p), rho_p)
          ALLOCATE(V_lr(SIZE(V_colin)))
          CALL compute_range_separation(rho_d%pot, rho_p%pot, V_d%nr1, V_d%nr2, V_d%nr3, &
               V_d%at, V_d%alat, coulomb_2d, alpha_gauss, V_colin, V_lr, SIZE(V_colin), &
               TRIM(edi_prefix))
          IF (ALLOCATED(rho_d%pot)) DEALLOCATE(rho_d%pot)
          IF (ALLOCATED(rho_p%pot)) DEALLOCATE(rho_p%pot)
          DEALLOCATE(V_lr)
        END BLOCK
     ENDIF
  ELSE
     IF (ionode) WRITE(stdout,'(/,5X,A)') &
          'potfile_d/potfile_p not set: skipping supercell-potential load (audit only).'
  ENDIF

  IF (ionode) THEN
     WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
     WRITE(stdout,'(5X,A)')   'EDT Phase 0 complete: inputs read, partition built. No solve yet.'
     WRITE(stdout,'(5X,A)')   REPEAT('=',64)
  ENDIF

  IF (ALLOCATED(irvec))     DEALLOCATE(irvec)
  IF (ALLOCATED(ndegen))    DEALLOCATE(ndegen)
  IF (ALLOCATED(chw))       DEALLOCATE(chw)
  IF (ALLOCATED(eig_interp))DEALLOCATE(eig_interp)
  CALL environment_end('EDT')
  CALL mp_global_end()
END PROGRAM edt
