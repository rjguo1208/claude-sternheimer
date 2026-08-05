MODULE edt_wannier
  !---------------------------------------------------------------------------
  !  Minimal, robust reader for EDI's `filukk` (combined Wannier rotation).
  !
  !  Why not EDI's read_filukk_edi?  That routine also reads the trailing
  !  lwindow / excluded_band(nbnd) blocks assuming the *current* nbnd.  For the
  !  T-matrix we deliberately use a high-band NSCF (nbnd=150) while filukk was
  !  written from the small Wannierization run (nbnd=17), so that block size
  !  mismatches and overruns the file.  Here we read only what we need —
  !  nbndep, nbndskip, ibndkept, and u_kc(nbndep,n_wannier,nkstot) — and stop.
  !  This sets u_mat / n_wannier / num_bands / nbndep / nbndskip / ibndkept in
  !  the shared EDI modules, exactly as the downstream interpolation expects.
  !---------------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: edt_read_filukk

CONTAINS

  SUBROUTINE edt_read_filukk(fname, nkstot, nks, nbndsub)
    USE io_global,   ONLY : ionode, ionode_id, stdout
    USE mp,          ONLY : mp_bcast
    USE mp_world,    ONLY : world_comm
    USE wann_common, ONLY : u_mat, u_mat_opt, n_wannier, num_bands, iknum, wann_centers
    USE global_var,  ONLY : nbndep, nbndskip, ibndkept
    USE parallelism, ONLY : fkbounds
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: nkstot, nks, nbndsub

    INTEGER :: iun, ik, ib, iw, nbe, nbs, lower, upper, ios
    COMPLEX(dp), ALLOCATABLE :: ukc(:,:,:)

    iun = 91
    iknum = nkstot
    n_wannier = nbndsub
    CALL mp_bcast(n_wannier, ionode_id, world_comm)

    IF (ionode) THEN
       OPEN(iun, FILE=TRIM(fname), FORM='formatted', STATUS='old', IOSTAT=ios)
       IF (ios /= 0) CALL errore('edt_read_filukk', 'cannot open '//TRIM(fname), 1)
       READ(iun, *) nbe, nbs
    ENDIF
    CALL mp_bcast(nbe, ionode_id, world_comm)
    CALL mp_bcast(nbs, ionode_id, world_comm)
    nbndep = nbe; nbndskip = nbs; num_bands = nbe

    IF (ALLOCATED(ibndkept)) DEALLOCATE(ibndkept)
    ALLOCATE(ibndkept(nbndep))
    IF (ionode) THEN
       DO ib = 1, nbndep
          READ(iun, *) ibndkept(ib)
       ENDDO
    ENDIF
    CALL mp_bcast(ibndkept, ionode_id, world_comm)

    ALLOCATE(ukc(nbndep, n_wannier, nkstot))
    ukc = (0.0_dp, 0.0_dp)
    IF (ionode) THEN
       DO ik = 1, nkstot
          DO ib = 1, nbndep
             DO iw = 1, n_wannier
                READ(iun, *, IOSTAT=ios) ukc(ib, iw, ik)
                IF (ios /= 0) CALL errore('edt_read_filukk', 'EOF/error reading u_kc', 1)
             ENDDO
          ENDDO
       ENDDO
       CLOSE(iun)        ! ignore trailing lwindow/excluded_band/centers (size depends on nbnd)
    ENDIF
    CALL mp_bcast(ukc, ionode_id, world_comm)

    CALL fkbounds(nkstot, lower, upper)
    IF (ALLOCATED(u_mat)) DEALLOCATE(u_mat)
    ALLOCATE(u_mat(nbndep, n_wannier, nks))
    u_mat(:, :, 1:nks) = ukc(:, :, lower:upper)
    IF (ALLOCATED(u_mat_opt)) DEALLOCATE(u_mat_opt)
    ALLOCATE(u_mat_opt(nbndep, nbndep, nks))
    u_mat_opt = (0.0_dp, 0.0_dp)
    IF (.NOT. ALLOCATED(wann_centers)) ALLOCATE(wann_centers(3, n_wannier))
    wann_centers = 0.0_dp
    DEALLOCATE(ukc)

    IF (ionode) WRITE(stdout, '(5X,A,I4,A,I4,A,I4,A,I5)') &
         'filukk: nbndep=', nbndep, '  nbndskip=', nbndskip, &
         '  n_wannier=', n_wannier, '  nkstot=', nkstot
  END SUBROUTINE edt_read_filukk

END MODULE edt_wannier
