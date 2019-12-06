module hydro_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream
  use amrex_constants_module, only: ZERO, HALF, ONE

  implicit none

contains

  CASTRO_FORT_DEVICE subroutine divu(lo, hi, &
                                     q, q_lo, q_hi, &
                                     dx, div, div_lo, div_hi) bind(C, name='divu')
    ! this computes the *node-centered* divergence

    use amrex_constants_module, only: FOURTH
    use castro_module, only: QU, QV, QW, QVAR

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    integer  :: i, j, k
    real(rt) :: ux, vy, wz, dxinv, dyinv, dzinv

    dxinv = ONE / dx(1)
    dyinv = ONE / dx(2)
    dzinv = ONE / dx(3)

    !$acc parallel loop gang vector collapse(3) deviceptr(div, q) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ux = FOURTH*( &
                  + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                  + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                  + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                  + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                  + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                  + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                  + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                  + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                  + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                  + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                  + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                  + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine divu



  CASTRO_FORT_DEVICE subroutine ctoprim(lo, hi, &
                                        u,     u_lo,   u_hi, &
                                        q,     q_lo,   q_hi, &
                                        qaux, qa_lo,  qa_hi) bind(c,name='ctoprim')

    use network, only: nspec, aion_inv, zion
    use eos_module, only: eos_t, eos_input_re, eos
    use castro_module, only: NVAR, URHO, UMX, UMZ, &
                             UEDEN, UEINT, UTEMP, &
                             QRHO, QU, QV, QW, UFS, &
                             QREINT, QPRES, QTEMP, QGAME, QFS, &
                             QVAR, QC, QGAMC, QDPDR, QDPDE, NQAUX, &
                             small_dens, dual_energy_eta1

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    integer  :: i, j, k, g
    integer  :: n, iq, ispec
    real(rt) :: kineng, rhoinv
    real(rt) :: vel(3)

    type (eos_t) :: eos_state

    !$acc parallel loop gang vector collapse(3) deviceptr(u, q, qaux) private(vel) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             q(i,j,k,QRHO) = u(i,j,k,URHO)
             rhoinv = ONE/q(i,j,k,QRHO)

             vel = u(i,j,k,UMX:UMZ) * rhoinv

             q(i,j,k,QU:QW) = vel

             ! Get the internal energy, which we'll use for
             ! determining the pressure.  We use a dual energy
             ! formalism. If (E - K) < eta1 and eta1 is suitably
             ! small, then we risk serious numerical truncation error
             ! in the internal energy.  Therefore we'll use the result
             ! of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (u(i,j,k,UEDEN) - kineng) / u(i,j,k,UEDEN) .gt. dual_energy_eta1) then
                q(i,j,k,QREINT) = (u(i,j,k,UEDEN) - kineng) * rhoinv
             else
                q(i,j,k,QREINT) = u(i,j,k,UEINT) * rhoinv
             endif

             q(i,j,k,QTEMP) = u(i,j,k,UTEMP)

             ! Load passively advected quatities into q
             do ispec = 1, nspec
                n  = UFS + ispec - 1
                iq = QFS + ispec - 1
                q(i,j,k,iq) = u(i,j,k,n) * rhoinv
             enddo

             ! get gamc, p, T, c, csml using q state
             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % abar = ONE / (sum(q(i,j,k,QFS:QFS+nspec-1) * aion_inv(:)))
             eos_state % zbar = eos_state % abar * (sum(q(i,j,k,QFS:QFS+nspec-1) * zion(:) * aion_inv(:)))

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e * q(i,j,k,QRHO)
             q(i,j,k,QPRES)  = eos_state % p
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE

             qaux(i,j,k,QDPDR)  = eos_state % dpdr_e
             qaux(i,j,k,QDPDE)  = eos_state % dpde
             qaux(i,j,k,QGAMC)  = eos_state % gam1
             qaux(i,j,k,QC   )  = eos_state % cs

          enddo
       enddo
    enddo

  end subroutine ctoprim



  CASTRO_FORT_DEVICE subroutine apply_av(lo, hi, idir, dx, &
                                         div, div_lo, div_hi, &
                                         u, u_lo, u_hi, &
                                         flux, f_lo, f_hi) bind(c, name="apply_av")

    use amrex_constants_module, only: FOURTH
    use castro_module, only: NVAR, UTEMP

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir

    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    integer :: i, j, k, n

    real(rt) :: div1

    real(rt), parameter :: difmag = 0.1d0

    !$acc parallel loop gang vector collapse(4) deviceptr(flux, u, div) async(acc_stream)
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if ( n == UTEMP ) cycle

                if (idir .eq. 1) then

                   div1 = FOURTH * (div(i,j,k  ) + div(i,j+1,k) + &
                                    div(i,j,k+1) + div(i,j+1,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (u(i,j,k,n) - u(i-1,j,k,n))

                else if (idir .eq. 2) then

                   div1 = FOURTH * (div(i,j,k  ) + div(i+1,j,k) + &
                                    div(i,j,k+1) + div(i+1,j,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (u(i,j,k,n) - u(i,j-1,k,n))

                else

                   div1 = FOURTH * (div(i,j  ,k) + div(i+1,j  ,k) + &
                                    div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (u(i,j,k,n) - u(i,j,k-1,n))

                end if

                flux(i,j,k,n) = flux(i,j,k,n) + dx(idir) * div1

             end do
          end do
       end do

    end do

  end subroutine apply_av

  
  CASTRO_FORT_DEVICE subroutine normalize_species_fluxes(lo, hi, flux, f_lo, f_hi) bind(c, name="normalize_species_fluxes")
    ! Normalize the fluxes of the mass fractions so that
    ! they sum to 0.  This is essentially the CMA procedure that is
    ! defined in Plewa & Muller, 1999, A&A, 342, 179.

    use network, only: nspec
    use amrex_fort_module, only: rt => amrex_real
    use castro_module, only: NVAR, URHO, UFS

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k, n
    real(rt) :: sum, fac

    !$acc parallel loop gang vector collapse(3) deviceptr(flux) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             sum = ZERO

             do n = UFS, UFS+nspec-1
                sum = sum + flux(i,j,k,n)
             end do

             if (sum .ne. ZERO) then
                fac = flux(i,j,k,URHO) / sum
             else
                fac = ONE
             end if

             do n = UFS, UFS+nspec-1
                flux(i,j,k,n) = flux(i,j,k,n) * fac
             end do

          end do
       end do
    end do

  end subroutine normalize_species_fluxes



  CASTRO_FORT_DEVICE subroutine store_flux(lo, hi, &
                                           flux_out, fo_lo, fo_hi, &
                                           flux_in, fi_lo, fi_hi, &
                                           area, a_lo, a_hi, &
                                           dt) bind(C, name="store_flux")

    use castro_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: fo_lo(3), fo_hi(3)
    integer,  intent(in   ) :: fi_lo(3), fi_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(inout) :: flux_out(fo_lo(1):fo_hi(1),fo_lo(2):fo_hi(2),fo_lo(3):fo_hi(3),NVAR)
    real(rt), intent(in   ) :: flux_in(fi_lo(1):fi_hi(1),fi_lo(2):fi_hi(2),fi_lo(3):fi_hi(3),NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    real(rt), intent(in), value :: dt

    integer :: i, j, k, n

    !$acc parallel loop gang vector collapse(4) deviceptr(flux_out, flux_in, area) async(acc_stream)
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux_out(i,j,k,n) = dt * flux_in(i,j,k,n) * area(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine store_flux

  

  CASTRO_FORT_DEVICE subroutine fill_hydro_source(lo, hi, &
                                                  u, u_lo, u_hi, &
                                                  q, q_lo, q_hi, &
                                                  source, sr_lo, sr_hi, &
                                                  flux1, flux1_lo, flux1_hi, &
                                                  flux2, flux2_lo, flux2_hi, &
                                                  flux3, flux3_lo, flux3_hi, &
                                                  qx, qx_lo, qx_hi, &
                                                  qy, qy_lo, qy_hi, &
                                                  qz, qz_lo, qz_hi, &
                                                  area1, area1_lo, area1_hi, &
                                                  area2, area2_lo, area2_hi, &
                                                  area3, area3_lo, area3_hi, &
                                                  vol, vol_lo, vol_hi, &
                                                  dx, dt) bind(C, name="fill_hydro_source")

    use castro_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, &
                             UEINT, UTEMP, NGDNV, QVAR, &
                             GDU, GDV, GDW, GDPRES

    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::     u_lo(3),     u_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) ::    sr_lo(3),    sr_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)

    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: source(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3),NVAR)
    real(rt), intent(in) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt), intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt), intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    integer :: i, j, g, k, n
    real(rt) :: volInv
    real(rt) :: pdivu

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update.

    !$acc parallel loop gang vector collapse(3) deviceptr(source, flux1, flux2, flux3, area1, area2, area3) &
    !$acc deviceptr(qx, qy, qz, vol) async(acc_stream)
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                source(i,j,k,n) = source(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) &
                     + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) &
                     + flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) &
                     ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   pdivu = HALF * (qx(i+1,j,k,GDPRES) + qx(i,j,k,GDPRES)) * &
                                  (qx(i+1,j,k,GDU) - qx(i,j,k,GDU)) / dx(1) + &
                           HALF * (qy(i,j+1,k,GDPRES) + qy(i,j,k,GDPRES)) * &
                                  (qy(i,j+1,k,GDV) - qy(i,j,k,GDV)) / dx(2) + &
                           HALF * (qz(i,j,k+1,GDPRES) + qz(i,j,k,GDPRES)) * &
                                  (qz(i,j,k+1,GDW) - qz(i,j,k,GDW)) / dx(3)

                   source(i,j,k,n) = source(i,j,k,n) - pdivu
                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine fill_hydro_source

end module hydro_module
