! Copyright (c) 2026, University of Michigan
!
! Kinematic transport tests for MPAS-Atmosphere: the dynamics are switched off and
! the transport mass fluxes are prescribed from an analytic flow, so that the
! production transport code can be run on test cases with known exact solutions.
!
module atm_kinematic

   use mpas_kind_types, only : RKIND, StrKIND
   use mpas_derived_types, only : block_type, mpas_pool_type, MPAS_LOG_ERR, MPAS_LOG_CRIT
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, &
                                  mpas_pool_get_config, mpas_pool_get_array
   use mpas_constants, only : pii
   use mpas_log, only : mpas_log_write

   implicit none

   private
   public :: atm_kinematic_check_config, atm_kinematic_set_fluxes


contains


   !-----------------------------------------------------------------------
   !  routine atm_kinematic_check_config
   !
   !> \brief Check the kinematic-test namelist settings at start-up
   !> \details Sets ierr = 1 if config_kinematic_case is unknown, or is
   !>  not 'off' while split transport is off.
   !-----------------------------------------------------------------------
   subroutine atm_kinematic_check_config(configs, ierr)

      implicit none

      type (mpas_pool_type), intent(in) :: configs
      integer, intent(inout) :: ierr

      character (len=StrKIND), pointer :: config_kinematic_case
      logical, pointer :: config_split_dynamics_transport

      call mpas_pool_get_config(configs, 'config_kinematic_case', config_kinematic_case)
      if (trim(config_kinematic_case) == 'off') return

      select case (trim(config_kinematic_case))
      case ('lauritzen_nondiv')
      case default
         call mpas_log_write('Unknown config_kinematic_case '''//trim(config_kinematic_case)//'''', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
      end select

      call mpas_pool_get_config(configs, 'config_split_dynamics_transport', config_split_dynamics_transport)
      if (.not. config_split_dynamics_transport) then
         call mpas_log_write('Kinematic transport requires config_split_dynamics_transport = true.', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
      end if

   end subroutine atm_kinematic_check_config


   !-----------------------------------------------------------------------
   !  routine atm_kinematic_set_fluxes
   !
   !> \brief Prescribe ruAvg, wwAvg and rho_zz for the current step
   !> \details Called once per step, before the transport loop of atm_srk3,
   !>  with t0 the time at the start of the step. Dispatches on
   !>  config_kinematic_case. With check = .true. (first step) the case
   !>  routine also checks its fluxes against the analytic wind.
   !-----------------------------------------------------------------------
   subroutine atm_kinematic_set_fluxes(block, t0, dt, check)

      implicit none

      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: t0, dt
      logical, intent(in) :: check

      character (len=StrKIND), pointer :: config_kinematic_case

      call mpas_pool_get_config(block % configs, 'config_kinematic_case', config_kinematic_case)

      select case (trim(config_kinematic_case))
      case ('lauritzen_nondiv')
         call kinematic_lauritzen_nondiv(block, t0, dt, check)
      case default
         call mpas_log_write('atm_kinematic_set_fluxes: unknown case '''//trim(config_kinematic_case)//'''', &
                             messageType=MPAS_LOG_CRIT)
      end select

   end subroutine atm_kinematic_set_fluxes


   !-----------------------------------------------------------------------
   !  routine kinematic_lauritzen_nondiv
   !
   !> \brief Fluxes for the non-divergent flow of Lauritzen et al. (2012, GMD 5, 887-901)
   !> \details
   !>  Deformational flow with a background solid-body rotation, eqs. (18)-(20),
   !>  with T = 12 days and k = 10a/T. The solid-body rotation term
   !>  (2 pi a/T) cos(lat) carries the pattern once around the sphere in T, so
   !>  the tracer crosses different cells and edges after the flow reverses at
   !>  T/2, and errors made in the first half are not undone by the reversal.
   !>
   !>  The streamfunction is averaged over [t0, t0+dt] with 3-point
   !>  Gauss-Legendre quadrature at the cell vertices, and the normal flux on
   !>  each edge is its difference along the edge, so the flux divergence of
   !>  every cell is zero to roundoff. rho_zz is set to 1 in both time levels
   !>  and wwAvg to 0, so every level is an independent copy of the 2D problem.
   !>  With check = .true. the fluxes are compared with the point-sampled
   !>  analytic normal wind at mid-step, and the run stops if they disagree
   !>  (for example, a sign error in the edge orientation).
   !-----------------------------------------------------------------------
   subroutine kinematic_lauritzen_nondiv(block, t0, dt, check)

      implicit none

      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: t0, dt
      logical, intent(in) :: check

      type (mpas_pool_type), pointer :: mesh, state, diag
      integer, pointer :: nVertLevels, nEdges, nVertices
      integer, dimension(:,:), pointer :: verticesOnEdge
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex, latEdge, lonEdge, angleEdge, dvEdge
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, wwAvg, rho_zz_1, rho_zz_2, u_2
      real (kind=RKIND), pointer :: sphere_radius

      real (kind=RKIND), parameter :: period = 12.0_RKIND * 86400.0_RKIND
      ! 3-point Gauss-Legendre rule on [-1,1]: nodes 0, +-sqrt(3/5); weights 8/9, 5/9 (halved here, so they sum to 1)
      real (kind=RKIND), dimension(3), parameter :: gl_x = (/ -0.7745966692414834_RKIND, 0.0_RKIND, &
                                                               0.7745966692414834_RKIND /)
      real (kind=RKIND), dimension(3), parameter :: gl_w = (/ 5.0_RKIND, 8.0_RKIND, 5.0_RKIND /) / 18.0_RKIND

      real (kind=RKIND), dimension(:), allocatable :: psi
      real (kind=RKIND) :: a, kdef, t, lamp, un, flux, err_max, un_max
      integer :: iEdge, iVtx, ig, v1, v2, n_fallback

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'latVertex', latVertex)
      call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_1, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_2, 2)
      call mpas_pool_get_array(state, 'u', u_2, 2)

      a = sphere_radius
      kdef = 10.0_RKIND * a / period

      ! Streamfunction (m^2/s) at the vertices, averaged over the step: u = -(1/a) dpsi/dlat, v = (1/(a cos(lat))) dpsi/dlon
      allocate(psi(nVertices+1))
      psi(:) = 0.0_RKIND
      do iVtx = 1, nVertices
         do ig = 1, 3
            t = t0 + 0.5_RKIND * dt * (1.0_RKIND + gl_x(ig))
            lamp = lonVertex(iVtx) - 2.0_RKIND * pii * t / period
            psi(iVtx) = psi(iVtx) + gl_w(ig) * a * ( kdef * sin(lamp)**2 * cos(latVertex(iVtx))**2 * cos(pii * t / period) &
                                                     - 2.0_RKIND * pii * a / period * sin(latVertex(iVtx)) )
         end do
      end do

      ! Normal flux on each edge. The tangent t = k x n points from verticesOnEdge(1) to verticesOnEdge(2),
      ! so u.n = (k x grad psi).n = -grad psi . t. Outer-halo edges with a missing vertex use the
      ! point-sampled wind; they are never used to update owned cells.
      n_fallback = 0
      err_max = 0.0_RKIND
      un_max = 0.0_RKIND
      do iEdge = 1, nEdges
         v1 = verticesOnEdge(1,iEdge)
         v2 = verticesOnEdge(2,iEdge)
         un = kinematic_normal_wind(latEdge(iEdge), lonEdge(iEdge), angleEdge(iEdge), t0 + 0.5_RKIND*dt, a, kdef, period)
         if (v1 >= 1 .and. v1 <= nVertices .and. v2 >= 1 .and. v2 <= nVertices) then
            flux = -(psi(v2) - psi(v1)) / dvEdge(iEdge)
            err_max = max(err_max, abs(flux - un))
            un_max = max(un_max, abs(un))
         else
            flux = un
            n_fallback = n_fallback + 1
         end if
         ruAvg(1:nVertLevels,iEdge) = flux
         u_2(1:nVertLevels,iEdge) = flux      ! diagnostic only: rho_zz = 1, so u = ruAvg (step-averaged wind)
      end do
      ruAvg(:,nEdges+1) = 0.0_RKIND

      wwAvg(:,:) = 0.0_RKIND
      rho_zz_1(:,:) = 1.0_RKIND
      rho_zz_2(:,:) = 1.0_RKIND

      if (check) then
         call mpas_log_write('Kinematic fluxes: max |u.n| = $r m/s, max |psi flux - point wind| = $r m/s, fallback edges = $i', &
                             realArgs=(/un_max, err_max/), intArgs=(/n_fallback/))
         if (err_max > 0.1_RKIND * un_max) then
            call mpas_log_write('Kinematic fluxes disagree with the analytic wind; check the edge orientation.', &
                                messageType=MPAS_LOG_CRIT)
         end if
      end if

      deallocate(psi)

   contains

      real (kind=RKIND) function kinematic_normal_wind(lat, lon, angle, t, a, kdef, period)
         real (kind=RKIND), intent(in) :: lat, lon, angle, t, a, kdef, period
         real (kind=RKIND) :: lamp, u, v
         lamp = lon - 2.0_RKIND * pii * t / period
         u = kdef * sin(lamp)**2 * sin(2.0_RKIND*lat) * cos(pii * t / period) + 2.0_RKIND * pii * a / period * cos(lat)
         v = kdef * sin(2.0_RKIND*lamp) * cos(lat) * cos(pii * t / period)
         kinematic_normal_wind = u * cos(angle) + v * sin(angle)
      end function kinematic_normal_wind

   end subroutine kinematic_lauritzen_nondiv

end module atm_kinematic
