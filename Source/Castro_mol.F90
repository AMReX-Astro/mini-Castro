! advection routines in support of method of lines integration

module mol_module

  implicit none

contains
  
#ifdef CUDA
  attributes(device) &
#endif
  subroutine mol_single_stage(time, &
                              lo, hi, domlo, domhi, &
                              uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                              uout, uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3, &
                              q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                              qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3, &
                              update, updt_l1, updt_l2, updt_l3, updt_h1, updt_h2, updt_h3, &
                              dx, dt, h, &
                              flux1, flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3, &
                              flux2, flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3, &
                              flux3, flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3, &
                              area1, area1_l1, area1_l2, area1_l3, area1_h1, area1_h2, area1_h3, &
                              area2, area2_l1, area2_l2, area2_l3, area2_h1, area2_h2, area2_h3, &
                              area3, area3_l1, area3_l2, area3_l3, area3_h1, area3_h2, area3_h3, &
                              vol, vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3, &
                              courno, verbose)

    use advection_util_module, only: compute_cfl, divu, normalize_species_fluxes, &
                                     ht
    use bl_constants_module, only: ZERO, HALF, ONE, FOURTH
    use flatten_module, only: uflaten
    use riemann_module, only: cmpflx
    use ppm_module, only: ppm_reconstruct
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                   UTEMP, UEINT, UMX, GDU, GDV, GDW, &
                                   QU, QV, QW, QPRES, NQAUX

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), verbose
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
    integer,  intent(in   ) :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    integer,  intent(in   ) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
    integer,  intent(in   ) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
    integer,  intent(in   ) :: updt_l1, updt_l2, updt_l3, updt_h1, updt_h2, updt_h3
    integer,  intent(in   ) :: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
    integer,  intent(in   ) :: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
    integer,  intent(in   ) :: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3
    integer,  intent(in   ) :: area1_l1, area1_l2, area1_l3, area1_h1, area1_h2, area1_h3
    integer,  intent(in   ) :: area2_l1, area2_l2, area2_l3, area2_h1, area2_h2, area2_h3
    integer,  intent(in   ) :: area3_l1, area3_l2, area3_l3, area3_h1, area3_h2, area3_h3
    integer,  intent(in   ) :: vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3

    real(rt), intent(in   ) :: uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3, NVAR)
    real(rt), intent(inout) :: uout(uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NVAR)
    real(rt), intent(inout) :: q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3, NQ)
    real(rt), intent(inout) :: qaux(qa_l1:qa_h1, qa_l2:qa_h2, qa_l3:qa_h3, NQAUX)
    real(rt), intent(inout) :: update(updt_l1:updt_h1, updt_l2:updt_h2, updt_l3:updt_h3, NVAR)
    real(rt), intent(inout) :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, flux1_l3:flux1_h3, NVAR)
    real(rt), intent(inout) :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, flux2_l3:flux2_h3, NVAR)
    real(rt), intent(inout) :: flux3(flux3_l1:flux3_h1, flux3_l2:flux3_h2, flux3_l3:flux3_h3, NVAR)
    real(rt), intent(in   ) :: area1(area1_l1:area1_h1, area1_l2:area1_h2, area1_l3:area1_h3)
    real(rt), intent(in   ) :: area2(area2_l1:area2_h1, area2_l2:area2_h2, area2_l3:area2_h3)
    real(rt), intent(in   ) :: area3(area3_l1:area3_h1, area3_l2:area3_h2, area3_l3:area3_h3)
    real(rt), intent(in   ) :: vol(vol_l1:vol_h1, vol_l2:vol_h2, vol_l3:vol_h3)
    real(rt), intent(in   ) :: dx(3), dt, time
    type(ht), intent(inout) :: h
    real(rt), intent(inout) :: courno

    integer :: ngf
    integer :: uin_lo(3), uin_hi(3)
    integer :: uout_lo(3), uout_hi(3)
    integer :: flux1_lo(3), flux1_hi(3)
    integer :: flux2_lo(3), flux2_hi(3)
    integer :: flux3_lo(3), flux3_hi(3)
    integer :: area1_lo(3), area1_hi(3)
    integer :: area2_lo(3), area2_hi(3)
    integer :: area3_lo(3), area3_hi(3)
    integer :: vol_lo(3), vol_hi(3)
    integer :: q_lo(3), q_hi(3)
    integer :: qa_lo(3), qa_hi(3)

    integer :: edge_lo(3), edge_hi(3)
    integer :: g_lo(3), g_hi(3)

    real(rt), parameter :: difmag = 0.1d0

    real(rt) :: div1
    integer :: i, j, k, n
    integer :: kc, km, kt, k3d

    ngf = 1

    q_lo = [ q_l1, q_l2, q_l3 ]
    q_hi = [ q_h1, q_h2, q_h3 ]

    qa_lo = [ qa_l1, qa_l2, qa_l3 ]
    qa_hi = [ qa_h1, qa_h2, qa_h3 ]

    uin_lo = [ uin_l1, uin_l2, uin_l3 ]
    uin_hi = [ uin_h1, uin_h2, uin_h3 ]

    uout_lo = [ uout_l1, uout_l2, uout_l3 ]
    uout_hi = [ uout_h1, uout_h2, uout_h3 ]

    flux1_lo = [ flux1_l1, flux1_l2, flux1_l3 ]
    flux1_hi = [ flux1_h1, flux1_h2, flux1_h3 ]

    flux2_lo = [ flux2_l1, flux2_l2, flux2_l3 ]
    flux2_hi = [ flux2_h1, flux2_h2, flux2_h3 ]

    flux3_lo = [ flux3_l1, flux3_l2, flux3_l3 ]
    flux3_hi = [ flux3_h1, flux3_h2, flux3_h3 ]

    area1_lo = [ area1_l1, area1_l2, area1_l3 ]
    area1_hi = [ area1_h1, area1_h2, area1_h3 ]

    area2_lo = [ area2_l1, area2_l2, area2_l3 ]
    area2_hi = [ area2_h1, area2_h2, area2_h3 ]

    area3_lo = [ area3_l1, area3_l2, area3_l3 ]
    area3_hi = [ area3_h1, area3_h2, area3_h3 ]

    vol_lo = [ vol_l1, vol_l2, vol_l3 ]
    vol_hi = [ vol_h1, vol_h2, vol_h3 ]

    g_lo = lo - ngf
    g_hi = hi + ngf

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             h%shk(i,j,k) = ZERO
          enddo
       enddo
    enddo

    ! Check if we have violated the CFL criterion.
    call compute_cfl(q, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi, &
                     lo, hi, dt, dx, courno)

    ! Compute flattening coefficient for slope calculations.

    call uflaten(g_lo, g_hi, q, q_lo, q_hi, h)

    do n = 1, NQ

       call ppm_reconstruct(q, q_lo, q_hi, h, lo, hi, dx, n)

       ! Construct the interface states -- this is essentially just a
       ! reshuffling of interface states from zone-center indexing to
       ! edge-centered indexing
       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                ! x-edges

                ! left state at i-1/2 interface
                h%qm(i,j,k,n,1) = h%sxp(i-1,j,k,n)

                ! right state at i-1/2 interface
                h%qp(i,j,k,n,1) = h%sxm(i,j,k,n)

                ! y-edges

                ! left state at j-1/2 interface
                h%qm(i,j,k,n,2) = h%syp(i,j-1,k,n)

                ! right state at j-1/2 interface
                h%qp(i,j,k,n,2) = h%sym(i,j,k,n)

                ! z-edges

                ! left state at k3d-1/2 interface
                h%qm(i,j,k,n,3) = h%szp(i,j,k-1,n)

                ! right state at k3d-1/2 interface
                h%qp(i,j,k,n,3) = h%szm(i,j,k,n)

             end do
          end do
       end do

    end do

    ! Compute F^x at kc (k3d)
    call cmpflx(flux1, flux1_lo, flux1_hi, &
                qaux, qa_lo, qa_hi, &
                h, 1, lo, [hi(1)+1, hi(2), hi(3)], domlo, domhi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             h%q1(i,j,k,:) = h%qint(i,j,k,:)
          end do
       end do
    end do

    ! Compute F^y at kc (k3d)
    call cmpflx(flux2, flux2_lo, flux2_hi, &
                qaux, qa_lo, qa_hi, &
                h, 2, lo, [hi(1), hi(2)+1, hi(3)], domlo, domhi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             h%q2(i,j,k,:) = h%qint(i,j,k,:)
          end do
       end do
    end do

    ! Compute F^z at kc (k3d)

    call cmpflx(flux3, flux3_lo, flux3_hi, &
                qaux, qa_lo, qa_hi, &
                h, 3, lo, [hi(1), hi(2), hi(3)+1], domlo, domhi)

    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             h%q3(i,j,k,:) = h%qint(i,j,k,:)
          end do
       end do
    end do

    ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
    edge_lo = lo
    edge_hi = hi + 1
    call divu(lo,hi,q,q_lo,q_hi,dx,h%div,edge_lo,edge_hi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             h%pdivu(i,j,k) = &
                  HALF*(h%q1(i+1,j,k,GDPRES) + h%q1(i,j,k,GDPRES)) * &
                       (h%q1(i+1,j,k,GDU) - h%q1(i,j,k,GDU))/dx(1) + &
                  HALF*(h%q2(i,j+1,k,GDPRES) + h%q2(i,j,k,GDPRES)) * &
                       (h%q2(i,j+1,k,GDV) - h%q2(i,j,k,GDV))/dx(2) + &
                  HALF*(h%q3(i,j,k+1,GDPRES) + h%q3(i,j,k,GDPRES)) * &
                       (h%q3(i,j,k+1,GDW) - h%q3(i,j,k,GDW))/dx(3)
          enddo
       enddo
    enddo

    do n = 1, NVAR

       if ( n == UTEMP ) then
          flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
          flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
          flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO

       else
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1
                   div1 = FOURTH*(h%div(i,j,k) + h%div(i,j+1,k) + &
                                  h%div(i,j,k+1) + h%div(i,j+1,k+1))
                   div1 = difmag*min(ZERO,div1)

                   flux1(i,j,k,n) = flux1(i,j,k,n) + &
                        dx(1) * div1 * (uin(i,j,k,n)-uin(i-1,j,k,n))
                enddo
             enddo
          enddo

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                   div1 = FOURTH*(h%div(i,j,k) + h%div(i+1,j,k) + &
                                  h%div(i,j,k+1) + h%div(i+1,j,k+1))
                   div1 = difmag*min(ZERO,div1)

                   flux2(i,j,k,n) = flux2(i,j,k,n) + &
                        dx(2) * div1 * (uin(i,j,k,n)-uin(i,j-1,k,n))
                enddo
             enddo
          enddo

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   div1 = FOURTH*(h%div(i,j,k) + h%div(i+1,j,k) + &
                                  h%div(i,j+1,k) + h%div(i+1,j+1,k))
                   div1 = difmag*min(ZERO,div1)

                   flux3(i,j,k,n) = flux3(i,j,k,n) + &
                        dx(3) * div1 * (uin(i,j,k,n)-uin(i,j,k-1,n))
                enddo
             enddo
          enddo

       endif

    enddo

    call normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
                                  flux2,flux2_lo,flux2_hi, &
                                  flux3,flux3_lo,flux3_hi, &
                                  lo,hi)

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                update(i,j,k,n) = update(i,j,k,n) + &
                     (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                      flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                      flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - h%pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do n = 1, NVAR
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux3(i,j,k,n) = dt * flux3(i,j,k,n) * area3(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine mol_single_stage

end module mol_module
