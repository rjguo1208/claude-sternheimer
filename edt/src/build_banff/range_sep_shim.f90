MODULE range_sep
  ! Banff-port shim: the v8d EDI tree dropped the range-separation module.
  ! EDT reaches this code only when the namelist switches range_sep on (it is
  ! off in every stage-A run); calling it aborts with a clear message.
  IMPLICIT NONE
CONTAINS
  SUBROUTINE compute_range_separation(rhod, rhop, nr1, nr2, nr3, at, alat, &
                                      c2d, agauss, vsr, vlr, n, prefix)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)         :: rhod(*), rhop(*), vsr(*), vlr(*)
    INTEGER,  INTENT(IN) :: nr1, nr2, nr3, n
    REAL(DP), INTENT(IN) :: at(3,3), alat, agauss
    LOGICAL,  INTENT(IN) :: c2d
    CHARACTER(LEN=*), INTENT(IN) :: prefix
    CALL errore('compute_range_separation', &
                'range separation is not available in the Banff v8d-linked build', 1)
  END SUBROUTINE compute_range_separation
END MODULE range_sep
