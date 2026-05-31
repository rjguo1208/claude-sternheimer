MODULE edt_sternheimer
  !---------------------------------------------------------------------------
  !  Stage C (plan.md §6): the bare per-k' Sternheimer solve.
  !
  !  Solve   [ (H0(k') - omega0) + alpha*P_act ] |chi> = |s>     (s in rest)
  !  by projected CG using QE's h_psi as the matvec.  For this MoS2 case the
  !  rest manifold lies entirely above omega0 (=VBM), so (H0-omega0) is
  !  positive-definite on rest and real CG converges; the alpha*P_act shift
  !  lifts the active eigenvalues so the solution stays in rest.
  !  Then  <s|chi> = Σ_{r∈R}|<r|s>|²/(ε_r-omega0) = -ΔṼ_chan, giving the rest
  !  dressing with ALL bands summed implicitly (the point vs the explicit sum).
  !
  !  Prerequisite: h_psi needs vrs (set here once) + per-k g2kin/vkb.  The
  !  eigenvalue gate <ψ_nk|H0|ψ_nk> == ε_nk validates the whole setup.
  !---------------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: edt_set_vrs, hpsi_setup_k, test_hpsi_eigen, solve_rest_cg, rest_channel_compare

CONTAINS

  SUBROUTINE edt_set_vrs()
    !! Build vrs = vltot + v_Hxc on the smooth grid for h_psi (post read_file).
    USE scf,       ONLY : vrs, vltot, v, kedtau
    USE lsda_mod,  ONLY : nspin
    USE gvecs,     ONLY : doublegrid
    USE fft_base,  ONLY : dffts
    USE io_global, ONLY : ionode, stdout
    IMPLICIT NONE
    CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dffts%nnr, nspin, doublegrid)
    IF (ionode) WRITE(stdout,'(5X,A)') 'h_psi potential vrs set (set_vrs).'
  END SUBROUTINE edt_set_vrs


  SUBROUTINE hpsi_setup_k(ik)
    !! Per-k setup so h_psi(...) evaluates H0(ik): kinetic g2kin + KB projectors vkb.
    USE wvfct,     ONLY : current_k
    USE klist,     ONLY : ngk, igk_k, xk
    USE uspp,      ONLY : vkb, nkb
    USE uspp_init, ONLY : init_us_2
    USE lsda_mod,  ONLY : current_spin, isk
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ik
    current_k = ik
    current_spin = isk(ik)
    CALL g2_kin(ik)
    IF (nkb > 0) CALL init_us_2(ngk(ik), igk_k(1,ik), xk(1,ik), vkb)
  END SUBROUTINE hpsi_setup_k


  SUBROUTINE test_hpsi_eigen(ik, blist, nb)
    !! Gate: <ψ_{b,ik}|H0|ψ_{b,ik}> must equal ε_{b,ik} (validates vrs/vkb/g2kin/h_psi).
    USE io_global,        ONLY : ionode, stdout
    USE wvfct,            ONLY : npwx, nbnd, et
    USE klist,            ONLY : ngk
    USE uspp,             ONLY : nkb
    USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type
    USE noncollin_module, ONLY : npol
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE io_files,         ONLY : restart_dir
    USE constants,        ONLY : rytoev
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ik, nb, blist(nb)
    COMPLEX(dp), ALLOCATABLE :: evc_k(:,:), psi(:,:), hpsi(:,:)
    REAL(dp) :: e_h, maxd
    INTEGER  :: i, npw, n2

    IF (.NOT. ionode) RETURN
    npw = ngk(ik)
    n2  = npwx*npol
    CALL hpsi_setup_k(ik)

    ALLOCATE(evc_k(n2,nbnd))
    CALL read_collected_wfc(restart_dir(), ik, evc_k)
    ALLOCATE(psi(n2,nb), hpsi(n2,nb))
    DO i = 1, nb
       psi(:,i) = evc_k(:, blist(i))
    ENDDO
    CALL allocate_bec_type(nkb, nb, becp)
    CALL h_psi(npwx, npw, nb, psi, hpsi)
    CALL deallocate_bec_type(becp)

    WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
    WRITE(stdout,'(5X,A,I4)') 'P2b gate — <psi|H0|psi> vs eigenvalue, ik=', ik
    WRITE(stdout,'(7X,A6,2A16,A14)') 'band','et (eV)','<H0> (eV)','diff (eV)'
    maxd = 0.0_dp
    DO i = 1, nb
       e_h = DBLE(SUM(CONJG(psi(1:npw,i))*hpsi(1:npw,i)))
       IF (npol == 2) e_h = e_h + DBLE(SUM(CONJG(psi(npwx+1:npwx+npw,i))*hpsi(npwx+1:npwx+npw,i)))
       maxd = MAX(maxd, ABS(e_h - et(blist(i),ik))*rytoev)
       WRITE(stdout,'(7X,I6,2F16.6,ES14.3)') blist(i), et(blist(i),ik)*rytoev, &
            e_h*rytoev, (e_h - et(blist(i),ik))*rytoev
    ENDDO
    WRITE(stdout,'(5X,A,ES12.4,A)') 'max |<H0> - et| = ', maxd, ' eV'
    IF (maxd < 1.0d-6) THEN
       WRITE(stdout,'(5X,A)') 'PASS: h_psi reproduces NSCF eigenvalues -> H0 setup correct.'
    ELSE
       WRITE(stdout,'(5X,A)') 'FAIL: h_psi != et (check vrs / vkb / g2kin setup).'
    ENDIF
    WRITE(stdout,'(5X,A)') REPEAT('=',64)
    FLUSH(stdout)

    DEALLOCATE(evc_k, psi, hpsi)
  END SUBROUTINE test_hpsi_eigen


  SUBROUTINE solve_rest_cg(kp, s, omega0_ry, alpha_ry, nact, evc_act, thr, maxiter, &
                            chi, iters, resid)
    !! Preconditioned CG for [ (H0-omega0) + alpha*P_act ] chi = s  (A pos-def on rest).
    !! Jacobi preconditioner 1/max(g2kin,1).  Scalar (npol=1).  h_psi must be set
    !! up for kp first (hpsi_setup_k).
    USE wvfct,            ONLY : npwx, g2kin
    USE klist,            ONLY : ngk
    USE uspp,             ONLY : nkb
    USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type
    USE noncollin_module, ONLY : npol
    IMPLICIT NONE
    INTEGER,  INTENT(IN)  :: kp, nact, maxiter
    COMPLEX(dp), INTENT(IN)  :: s(:), evc_act(:,:)
    REAL(dp), INTENT(IN)  :: omega0_ry, alpha_ry, thr
    COMPLEX(dp), INTENT(OUT) :: chi(:)
    INTEGER,  INTENT(OUT) :: iters
    REAL(dp), INTENT(OUT) :: resid
    COMPLEX(dp), ALLOCATABLE :: r(:), p(:), Ap(:), z(:), hp(:,:), ptmp(:,:)
    REAL(dp), ALLOCATABLE :: prec(:)
    REAL(dp) :: rz, rz_new, pAp, ac, bc
    COMPLEX(dp) :: ovl
    INTEGER :: npw, it, a, ig, n2

    npw = ngk(kp); n2 = npwx*npol
    ALLOCATE(r(n2), p(n2), Ap(n2), z(n2), hp(n2,1), ptmp(n2,1), prec(npw))
    DO ig = 1, npw
       prec(ig) = 1.0_dp / MAX(g2kin(ig), 1.0_dp)
    ENDDO
    CALL allocate_bec_type(nkb, 1, becp)

    chi = (0.0_dp,0.0_dp)
    r = (0.0_dp,0.0_dp); r(1:npw) = s(1:npw)
    z = (0.0_dp,0.0_dp); z(1:npw) = prec(1:npw) * r(1:npw)
    p = z
    rz = DBLE(SUM(CONJG(r(1:npw))*z(1:npw)))
    iters = 0; resid = SQRT(DBLE(SUM(CONJG(r(1:npw))*r(1:npw))))

    DO it = 1, maxiter
       ! Ap = (H0-omega0) p + alpha P_act p
       ptmp(:,1) = p
       CALL h_psi(npwx, npw, 1, ptmp, hp)
       Ap = (0.0_dp,0.0_dp)
       Ap(1:npw) = hp(1:npw,1) - omega0_ry*p(1:npw)
       DO a = 1, nact
          ovl = SUM(CONJG(evc_act(1:npw,a))*p(1:npw))
          Ap(1:npw) = Ap(1:npw) + alpha_ry*ovl*evc_act(1:npw,a)
       ENDDO
       pAp = DBLE(SUM(CONJG(p(1:npw))*Ap(1:npw)))
       ac = rz / pAp
       chi(1:npw) = chi(1:npw) + ac*p(1:npw)
       r(1:npw)   = r(1:npw)   - ac*Ap(1:npw)
       resid = SQRT(DBLE(SUM(CONJG(r(1:npw))*r(1:npw))))
       iters = it
       IF (resid < thr) EXIT
       z(1:npw) = prec(1:npw) * r(1:npw)
       rz_new = DBLE(SUM(CONJG(r(1:npw))*z(1:npw)))
       bc = rz_new / rz
       p(1:npw) = z(1:npw) + bc*p(1:npw)
       rz = rz_new
    ENDDO

    CALL deallocate_bec_type(becp)
    DEALLOCATE(r, p, Ap, z, hp, ptmp, prec)
  END SUBROUTINE solve_rest_cg


  SUBROUTINE rest_channel_compare(ki, isrc, kp, q_cryst, omega0_ry, win_min_ry, win_max_ry, &
                                   nbndskip_in, thr, cutoffs, ncut)
    !! At channel kp: build s=QΔV|isrc,ki>, then compare the rest dressing from
    !! (a) the explicit rest-band sum (cutoff convergence) and (b) the Sternheimer
    !! CG solve <s|χ> (all bands implicit).  The Sternheimer value is the
    !! converged limit the explicit sum approaches — the T2 validation.
    USE io_global,        ONLY : ionode, stdout
    USE wvfct,            ONLY : npwx, nbnd, et
    USE klist,            ONLY : ngk, igk_k, xk
    USE fft_base,         ONLY : dffts
    USE fft_interfaces,   ONLY : invfft
    USE noncollin_module, ONLY : npol
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE io_files,         ONLY : restart_dir
    USE constants,        ONLY : rytoev
    USE edic_mod,         ONLY : V_d, V_p
    USE edt_source,       ONLY : build_source_ket, count_nkb, make_coeff
    USE edt_partition,    ONLY : apply_Qproj
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: ki, isrc, kp, nbndskip_in, ncut, cutoffs(ncut)
    REAL(dp), INTENT(IN) :: q_cryst(3), omega0_ry, win_min_ry, win_max_ry, thr

    COMPLEX(dp), ALLOCATABLE :: evc_ki(:,:), evc_kp(:,:), u_src(:), psic(:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d_ki(:,:), vkb_p_ki(:,:), bec_d(:), bec_p(:)
    COMPLEX(dp), ALLOCATABLE :: coeff_d(:), coeff_p(:), zeta(:), chi(:), evc_act(:,:)
    COMPLEX(dp), ALLOCATABLE :: dv(:)
    COMPLEX(dp) :: ovlp
    REAL(dp) :: eR, alpha_ry, resid, dv_stern
    INTEGER :: nkb_d, nkb_p, ig, r, c, npw_ki, npw_kp, nact, iters, n2

    IF (.NOT. ionode) RETURN
    npw_ki = ngk(ki); npw_kp = ngk(kp); n2 = npwx*npol

    ALLOCATE(evc_ki(n2,nbnd), evc_kp(n2,nbnd))
    CALL read_collected_wfc(restart_dir(), ki, evc_ki)
    CALL read_collected_wfc(restart_dir(), kp, evc_kp)

    ALLOCATE(u_src(dffts%nnr), psic(dffts%nnr))
    psic = (0.0_dp,0.0_dp)
    DO ig = 1, npw_ki
       psic(dffts%nl(igk_k(ig,ki))) = evc_ki(ig, isrc)
    ENDDO
    CALL invfft('Wave', psic, dffts); u_src = psic

    CALL count_nkb(V_d%nat, V_d%ityp, V_d%ntyp, nkb_d)
    CALL count_nkb(V_p%nat, V_p%ityp, V_p%ntyp, nkb_p)
    ALLOCATE(vkb_d_ki(npwx,nkb_d), vkb_p_ki(npwx,nkb_p), bec_d(nkb_d), bec_p(nkb_p), &
             coeff_d(nkb_d), coeff_p(nkb_p))
    CALL get_betavkb(npw_ki, igk_k(1,ki), xk(1,ki), vkb_d_ki, V_d%nat, V_d%ityp, V_d%tau, nkb_d)
    CALL get_betavkb(npw_ki, igk_k(1,ki), xk(1,ki), vkb_p_ki, V_p%nat, V_p%ityp, V_p%tau, nkb_p)
    DO c = 1, nkb_d
       bec_d(c) = SUM(CONJG(vkb_d_ki(1:npw_ki,c))*evc_ki(1:npw_ki,isrc))
    ENDDO
    DO c = 1, nkb_p
       bec_p(c) = SUM(CONJG(vkb_p_ki(1:npw_ki,c))*evc_ki(1:npw_ki,isrc))
    ENDDO
    CALL make_coeff(V_d%nat, V_d%ityp, V_d%ntyp, nkb_d, bec_d, coeff_d)
    CALL make_coeff(V_p%nat, V_p%ityp, V_p%ntyp, nkb_p, bec_p, coeff_p)

    ! source ket on channel kp
    ALLOCATE(zeta(npwx))
    CALL build_source_ket(kp, q_cryst, u_src, coeff_d, nkb_d, coeff_p, nkb_p, zeta)

    ! active states at kp + project source onto rest
    nact = 0
    DO r = nbndskip_in+1, nbnd
       eR = et(r,kp)
       IF (eR >= win_min_ry .AND. eR <= win_max_ry) nact = nact + 1
    ENDDO
    ALLOCATE(evc_act(n2,MAX(nact,1)))
    c = 0
    DO r = nbndskip_in+1, nbnd
       eR = et(r,kp)
       IF (eR >= win_min_ry .AND. eR <= win_max_ry) THEN
          c = c + 1; evc_act(:,c) = evc_kp(:,r)
       ENDIF
    ENDDO
    CALL apply_Qproj(npw_kp, npol, nact, evc_act, zeta)

    ! (a) explicit rest-band spectral sum with cutoff convergence
    ALLOCATE(dv(ncut)); dv = (0.0_dp,0.0_dp)
    DO r = nbndskip_in+1, nbnd
       eR = et(r,kp)
       IF (eR >= win_min_ry .AND. eR <= win_max_ry) CYCLE
       ovlp = SUM(CONJG(evc_kp(1:npw_kp,r))*zeta(1:npw_kp))
       DO c = 1, ncut
          IF (r <= cutoffs(c)) dv(c) = dv(c) + ABS(ovlp)**2/(omega0_ry - eR)
       ENDDO
    ENDDO

    ! (b) Sternheimer CG solve (all bands implicit)
    CALL hpsi_setup_k(kp)
    alpha_ry = 2.0_dp*(omega0_ry - win_min_ry)
    ALLOCATE(chi(npwx))
    CALL solve_rest_cg(kp, zeta, omega0_ry, alpha_ry, nact, evc_act, thr, 500, chi, iters, resid)
    ! <s|chi> = Σ_r |<r|s>|²/(ε_r-ω0) = -ΔṼ_chan
    dv_stern = -DBLE(SUM(CONJG(zeta(1:npw_kp))*chi(1:npw_kp)))

    WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
    WRITE(stdout,'(5X,A)') 'P2b — Sternheimer vs explicit rest dressing (T2), one channel'
    WRITE(stdout,'(5X,A,I4,A,I4,A,I4,A,3F8.4)') 'isrc=',isrc,' ki=',ki,' kp=',kp,' q=',q_cryst
    WRITE(stdout,'(5X,A,I4,A,ES10.2,A,I5)') 'CG iters=',iters,'  residual=',resid,'  nact=',nact
    WRITE(stdout,'(5X,A)') '  explicit cutoff   ΔṼ_chan (Ry)'
    DO c = 1, ncut
       WRITE(stdout,'(9X,I4,5X,ES16.8)') cutoffs(c), DBLE(dv(c))
    ENDDO
    WRITE(stdout,'(5X,A,ES16.8,A)') 'Sternheimer (all bands) = ', dv_stern, ' Ry'
    WRITE(stdout,'(5X,A,ES12.4)')   'Sternheimer - explicit(full) = ', dv_stern - DBLE(dv(ncut))
    WRITE(stdout,'(5X,A)') '  (explicit should approach Sternheimer; residual = high-band tail beyond cutoff)'
    WRITE(stdout,'(5X,A)') REPEAT('=',64)
    FLUSH(stdout)

    DEALLOCATE(evc_ki,evc_kp,u_src,psic,vkb_d_ki,vkb_p_ki,bec_d,bec_p,coeff_d,coeff_p)
    DEALLOCATE(zeta,chi,evc_act,dv)
  END SUBROUTINE rest_channel_compare

END MODULE edt_sternheimer
