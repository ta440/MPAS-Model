! Copyright (c) 2026, University of Michigan
!
! Kinematic transport tests for MPAS-Atmosphere: the dynamics are switched off and
! the transport mass fluxes are prescribed from an analytic flow, so that the
! production transport code can be run on test cases with known exact solutions.
!
module atm_kinematic

   use mpas_kind_types, only : RKIND, StrKIND
   use mpas_derived_types, only : domain_type, block_type, mpas_pool_type, MPAS_LOG_ERR, MPAS_LOG_CRIT
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, &
                                  mpas_pool_get_config, mpas_pool_get_array
   use mpas_constants, only : pii, rgas
   use mpas_log, only : mpas_log_write

   implicit none

   ! Provides definition of halo_exchange_routine
#include "mpas_halo_interface.inc"

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
      case ('lauritzen_nondiv', 'lauritzen_div', 'dcmip2012_1_2')
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
   !>  routine also checks its fluxes against the analytic wind. Cases in which
   !>  the density changes exchange the halo of rho_zz (time level 2) through
   !>  the halo group 'dynamics:rho_zz'.
   !-----------------------------------------------------------------------
   subroutine atm_kinematic_set_fluxes(domain, block, t0, dt, check, exchange_halo_group, config_gpu_aware_mpi)

      implicit none

      type (domain_type), intent(inout) :: domain
      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: t0, dt
      logical, intent(in) :: check
      procedure (halo_exchange_routine) :: exchange_halo_group
      logical, intent(in) :: config_gpu_aware_mpi

      character (len=StrKIND), pointer :: config_kinematic_case

      call mpas_pool_get_config(block % configs, 'config_kinematic_case', config_kinematic_case)

      select case (trim(config_kinematic_case))
      case ('lauritzen_nondiv')
         call kinematic_lauritzen_nondiv(block, t0, dt, check)
      case ('lauritzen_div')
         call kinematic_lauritzen_div(domain, block, t0, dt, check, exchange_halo_group, config_gpu_aware_mpi)
      case ('dcmip2012_1_2')
         call kinematic_dcmip2012_1_2(block, t0, dt, check)
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


   !-----------------------------------------------------------------------
   !  routine kinematic_lauritzen_div
   !
   !> \brief Fluxes for the divergent flow of Lauritzen et al. (2012, GMD 5, 887-901)
   !> \details
   !>  Deformational divergent flow with a background solid-body rotation,
   !>  eqs. (21)-(22), T = 12 days. In the frame moving with the rotation the
   !>  velocity is a fixed field times cos(pi t/T), so parcels and the density
   !>  return to their initial state at T.
   !>
   !>  The normal velocity on each edge is averaged over the edge (a great-circle
   !>  arc, on which the unit normal is constant) and over [t0, t0+dt] with 3x3
   !>  Gauss-Legendre points. The density is advanced with the stages of the
   !>  MPAS RK3 (dt/2, dt/2, dt, or dt/3, dt/2, dt for config_time_integration_order
   !>  = 3), as in the dynamics: ruAvg is the mass flux of the final stage, with
   !>  the edge density from the second stage, and rho_zz at the new time level is
   !>  updated from ruAvg with the same arithmetic as atm_advance_scalars_mono_work,
   !>  so a uniform tracer stays uniform. wwAvg = 0, so every level is a copy of
   !>  the 2D problem. On the first step (check = .true.) rho_zz is set to 1 and the
   !>  fluxes are compared with the point-sampled normal wind at the edge midpoints.
   !-----------------------------------------------------------------------
   subroutine kinematic_lauritzen_div(domain, block, t0, dt, check, exchange_halo_group, config_gpu_aware_mpi)

      implicit none

      type (domain_type), intent(inout) :: domain
      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: t0, dt
      logical, intent(in) :: check
      procedure (halo_exchange_routine) :: exchange_halo_group
      logical, intent(in) :: config_gpu_aware_mpi

      type (mpas_pool_type), pointer :: mesh, state, diag
      integer, pointer :: nVertLevels, nEdges, nCells, nVertices
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnEdge, cellsOnEdge, edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, latEdge, lonEdge, angleEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge, invAreaCell, rdzw
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, wwAvg, rho_zz_1, rho_zz_2, u_2
      real (kind=RKIND), pointer :: sphere_radius
      integer, pointer :: config_time_integration_order

      real (kind=RKIND), parameter :: period = 12.0_RKIND * 86400.0_RKIND
      ! 3-point Gauss-Legendre rule on [-1,1]: nodes 0, +-sqrt(3/5); weights 8/9, 5/9 (halved here, so they sum to 1)
      real (kind=RKIND), dimension(3), parameter :: gl_x = (/ -0.7745966692414834_RKIND, 0.0_RKIND, &
                                                               0.7745966692414834_RKIND /)
      real (kind=RKIND), dimension(3), parameter :: gl_w = (/ 5.0_RKIND, 8.0_RKIND, 5.0_RKIND /) / 18.0_RKIND

      real (kind=RKIND), dimension(:), allocatable :: un_avg           ! (nEdges+1) edge- and step-mean u.n
      real (kind=RKIND), dimension(3) :: stage_frac
      logical, dimension(:), allocatable :: cell_ok                    ! all edges of the cell are in the block
      real (kind=RKIND), dimension(3) :: p1, p2, pq, nrm
      real (kind=RKIND) :: a, t, s, arc, un, un_mid, err_max, un_max
      integer :: iCell, iEdge, i, k, ig, is, v1, v2, c1, c2, n_fallback

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_1, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_2, 2)
      call mpas_pool_get_array(state, 'u', u_2, 2)

      call mpas_pool_get_config(block % configs, 'config_time_integration_order', config_time_integration_order)

      a = sphere_radius

      ! Stage fractions of the RK3 in atm_srk3
      stage_frac(:) = (/ 0.5_RKIND, 0.5_RKIND, 1.0_RKIND /)
      if (config_time_integration_order == 3) stage_frac(1) = 1.0_RKIND / 3.0_RKIND

      ! 2D cases start from rho_zz = 1; afterwards the time-level shift carries rho_zz forward
      if (check) rho_zz_1(:,:) = 1.0_RKIND

      !
      ! Edge- and step-mean normal velocity
      !
      allocate(un_avg(nEdges+1))
      un_avg(:) = 0.0_RKIND
      n_fallback = 0
      err_max = 0.0_RKIND
      un_max = 0.0_RKIND
      do iEdge = 1, nEdges
         v1 = verticesOnEdge(1,iEdge)
         v2 = verticesOnEdge(2,iEdge)
         un_mid = point_normal_wind(latEdge(iEdge), lonEdge(iEdge), angleEdge(iEdge), t0 + 0.5_RKIND*dt)
         if (v1 >= 1 .and. v1 <= nVertices .and. v2 >= 1 .and. v2 <= nVertices) then
            p1 = (/ xVertex(v1), yVertex(v1), zVertex(v1) /)
            p2 = (/ xVertex(v2), yVertex(v2), zVertex(v2) /)
            p1 = p1 / norm2(p1)
            p2 = p2 / norm2(p2)
            nrm = (/ p2(2)*p1(3) - p2(3)*p1(2), p2(3)*p1(1) - p2(1)*p1(3), p2(1)*p1(2) - p2(2)*p1(1) /)
            nrm = nrm / norm2(nrm)                       ! n = (v2 x v1)/|v2 x v1|, constant along the arc
            arc = acos(max(-1.0_RKIND, min(1.0_RKIND, dot_product(p1, p2))))
            un = 0.0_RKIND
            do ig = 1, 3
               t = t0 + 0.5_RKIND * dt * (1.0_RKIND + gl_x(ig))
               do is = 1, 3
                  s = 0.5_RKIND * (1.0_RKIND + gl_x(is))
                  pq = (sin((1.0_RKIND - s) * arc) * p1 + sin(s * arc) * p2) / sin(arc)
                  un = un + gl_w(ig) * gl_w(is) * cartesian_normal_wind(pq, nrm, t)
               end do
            end do
            err_max = max(err_max, abs(un - un_mid))
            un_max = max(un_max, abs(un_mid))
         else
            un = un_mid                                  ! outer-halo edge; never updates an owned cell
            n_fallback = n_fallback + 1
         end if
         un_avg(iEdge) = un
         u_2(1:nVertLevels,iEdge) = un                   ! diagnostic only: step-mean wind
      end do

      allocate(cell_ok(nCells+1))
      cell_ok(:) = .false.
      do iCell = 1, nCells
         cell_ok(iCell) = all(edgesOnCell(1:nEdgesOnCell(iCell),iCell) <= nEdges)
      end do

      !
      ! Density through the three RK stages; ruAvg ends as the final-stage mass flux
      !
      wwAvg(:,:) = 0.0_RKIND
      ruAvg(:,nEdges+1) = 0.0_RKIND
      ! rho_zz(2) holds the stage densities; its halo is exchanged after every stage. The halo of
      ! rho_zz(1) is valid: it is 1 on the first step, and the exchanged rho_zz(2) of the last step after.
      rho_zz_2(:,:) = rho_zz_1(:,:)
      do ig = 1, 3
         do iEdge = 1, nEdges
            c1 = cellsOnEdge(1,iEdge)
            c2 = cellsOnEdge(2,iEdge)
            ruAvg(:,iEdge) = edge_mean(rho_zz_2, c1, c2) * un_avg(iEdge)    ! from rho^n, rho*, rho**
         end do
         if (ig < 3) then
            call update_density(rho_zz_1, rho_zz_2, stage_frac(ig) * dt)    ! rho*, then rho**
         else
            call update_density(rho_zz_1, rho_zz_2, dt)                     ! rho at t0 + dt
         end if
         call exchange_halo_group(domain, 'dynamics:rho_zz', config_gpu_aware_mpi)
      end do

      if (check) then
         call mpas_log_write('Kinematic fluxes (divergent): max |u.n| = $r m/s, max |edge mean - midpoint wind| = $r m/s,' &
                             //' fallback edges = $i', realArgs=(/un_max, err_max/), intArgs=(/n_fallback/))
         if (err_max > 0.1_RKIND * un_max) then
            call mpas_log_write('Kinematic fluxes disagree with the analytic wind; check the edge orientation.', &
                                messageType=MPAS_LOG_CRIT)
         end if
      end if

      deallocate(un_avg, cell_ok)

   contains

      ! Mean of rho over the two cells of an edge, on every level; one-sided if a cell is outside the block
      function edge_mean(rho, c1, c2) result(rho_e)
         real (kind=RKIND), dimension(:,:), intent(in) :: rho
         integer, intent(in) :: c1, c2
         real (kind=RKIND), dimension(nVertLevels) :: rho_e
         if (c1 <= nCells .and. c2 <= nCells) then
            rho_e(:) = 0.5_RKIND * (rho(1:nVertLevels,c1) + rho(1:nVertLevels,c2))
         else if (c1 <= nCells) then
            rho_e(:) = rho(1:nVertLevels,c1)
         else
            rho_e(:) = rho(1:nVertLevels,c2)
         end if
      end function edge_mean

      ! rho_new = rho_old - dt_stage * div(ruAvg), with the arithmetic of atm_advance_scalars_mono_work;
      ! cells with an edge outside the block keep rho_old
      subroutine update_density(rho_old, rho_new, dt_stage)
         real (kind=RKIND), dimension(:,:), intent(in) :: rho_old
         real (kind=RKIND), dimension(:,:), intent(inout) :: rho_new
         real (kind=RKIND), intent(in) :: dt_stage
         real (kind=RKIND), dimension(nVertLevels) :: div
         integer :: jCell, j, jEdge
         do jCell = 1, nCells
            if (.not. cell_ok(jCell)) then
               rho_new(1:nVertLevels,jCell) = rho_old(1:nVertLevels,jCell)
               cycle
            end if
            div(:) = 0.0_RKIND
            do j = 1, nEdgesOnCell(jCell)
               jEdge = edgesOnCell(j,jCell)
               div(:) = div(:) - edgesOnCell_sign(j,jCell) * ruAvg(1:nVertLevels,jEdge) * dvEdge(jEdge) * invAreaCell(jCell)
            end do
            rho_new(1:nVertLevels,jCell) = rho_old(1:nVertLevels,jCell) &
                 + dt_stage * (div(:) - rdzw(1:nVertLevels) * (wwAvg(2:nVertLevels+1,jCell) - wwAvg(1:nVertLevels,jCell)))
         end do
      end subroutine update_density

      ! Analytic wind (eqs. 21-22) at (lat, lon) and time t
      subroutine div_wind(lat, lon, t, u, v)
         real (kind=RKIND), intent(in) :: lat, lon, t
         real (kind=RKIND), intent(out) :: u, v
         real (kind=RKIND) :: lamp
         lamp = lon - 2.0_RKIND * pii * t / period
         u = -5.0_RKIND * a / period * sin(0.5_RKIND*lamp)**2 * sin(2.0_RKIND*lat) * cos(lat)**2 * cos(pii*t/period) &
             + 2.0_RKIND * pii * a / period * cos(lat)
         v = 2.5_RKIND * a / period * sin(lamp) * cos(lat)**3 * cos(pii*t/period)
      end subroutine div_wind

      ! u.n at the unit vector p for a Cartesian unit normal nrm
      real (kind=RKIND) function cartesian_normal_wind(p, nrm, t)
         real (kind=RKIND), dimension(3), intent(in) :: p, nrm
         real (kind=RKIND), intent(in) :: t
         real (kind=RKIND) :: lat, lon, u, v
         lat = asin(max(-1.0_RKIND, min(1.0_RKIND, p(3))))
         lon = atan2(p(2), p(1))
         call div_wind(lat, lon, t, u, v)
         cartesian_normal_wind = u * (-sin(lon) * nrm(1) + cos(lon) * nrm(2)) &
                               + v * (-sin(lat)*cos(lon) * nrm(1) - sin(lat)*sin(lon) * nrm(2) + cos(lat) * nrm(3))
      end function cartesian_normal_wind

      ! u.n at an edge midpoint, from angleEdge (the MPAS convention)
      real (kind=RKIND) function point_normal_wind(lat, lon, angle, t)
         real (kind=RKIND), intent(in) :: lat, lon, angle, t
         real (kind=RKIND) :: u, v
         call div_wind(lat, lon, t, u, v)
         point_normal_wind = u * cos(angle) + v * sin(angle)
      end function point_normal_wind

   end subroutine kinematic_lauritzen_div


   !-----------------------------------------------------------------------
   !  routine kinematic_dcmip2012_1_2
   !
   !> \brief Fluxes for DCMIP-2012 test 1-2, the Hadley-like meridional circulation
   !> \details
   !>  Eqs. (36)-(38) of the DCMIP-2012 test case document (v1.7); Kent et al.
   !>  (2014, QJRMS 140, 1279-1293), Sect. 4. Solid-body zonal flow plus K = 5
   !>  overturning cells in latitude-height that reverse at tau/2, tau = 1 day.
   !>  The mass flux is non-divergent in 3D and rho_zz is held fixed at its
   !>  initial (isothermal, hydrostatic) profile.
   !>
   !>  Zonal part: rho_zz times the vertex-streamfunction flux of u0 cos(lat),
   !>  divergence-free to roundoff. Meridional part: separable as
   !>  A(t) g(lat) cos(pi z/ztop); the step mean of cos(pi t/tau) and the layer
   !>  mean of cos(pi z/ztop) are exact, and the edge mean of g(lat) j.n uses
   !>  3-point Gauss-Legendre along the arc. wwAvg is then derived from the
   !>  discrete horizontal divergence of g j, so that the 3D mass divergence of
   !>  every cell is zero to roundoff and a uniform tracer stays uniform. With
   !>  check = .true. the fluxes are compared with the point-sampled analytic
   !>  winds, and the mass divergence of the owned cells is checked.
   !-----------------------------------------------------------------------
   subroutine kinematic_dcmip2012_1_2(block, t0, dt, check)

      implicit none

      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: t0, dt
      logical, intent(in) :: check

      type (mpas_pool_type), pointer :: mesh, state, diag
      integer, pointer :: nVertLevels, nEdges, nCells, nCellsSolve, nVertices
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnEdge, edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, latVertex, latEdge, angleEdge, latCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, invAreaCell, rdzw
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign, zgrid
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, wwAvg, rho_zz_1, rho_zz_2, u_2
      real (kind=RKIND), pointer :: sphere_radius

      real (kind=RKIND), parameter :: tau = 86400.0_RKIND, u0 = 40.0_RKIND, w0 = 0.15_RKIND
      real (kind=RKIND), parameter :: ztop = 12000.0_RKIND, p0 = 1.0e+05_RKIND, t_iso = 300.0_RKIND
      real (kind=RKIND), parameter :: kc = 5.0_RKIND                   ! K, the number of overturning cells
      ! 3-point Gauss-Legendre rule on [-1,1]: nodes 0, +-sqrt(3/5); weights 8/9, 5/9 (halved here, so they sum to 1)
      real (kind=RKIND), dimension(3), parameter :: gl_x = (/ -0.7745966692414834_RKIND, 0.0_RKIND, &
                                                               0.7745966692414834_RKIND /)
      real (kind=RKIND), dimension(3), parameter :: gl_w = (/ 5.0_RKIND, 8.0_RKIND, 5.0_RKIND /) / 18.0_RKIND

      real (kind=RKIND), dimension(:), allocatable :: psi, m_edge, vmer, s_int, rho_col
      real (kind=RKIND), dimension(3) :: p1, p2, pq, nrm
      real (kind=RKIND) :: a, rho0, amp, uz, mer, uz_mid, mer_mid, s, arc, lat, lon, cdiv, cdiv_exact
      real (kind=RKIND) :: err_u, max_u, err_w, max_w, res, res_scale, hdiv
      integer :: iCell, iEdge, iVtx, j, k, is, v1, v2, n_fallback

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'latVertex', latVertex)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_1, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_2, 2)
      call mpas_pool_get_array(state, 'u', u_2, 2)

      a = sphere_radius
      rho0 = p0 / (rgas * t_iso)

      if (check) then
         if (abs(zgrid(nVertLevels+1,1) - ztop) > 1.0e-6_RKIND * ztop) then
            call mpas_log_write('dcmip2012_1_2 needs a model top of 12000 m (config_ztop in the init).', &
                                messageType=MPAS_LOG_CRIT)
         end if
         do k = 1, nVertLevels
            if (maxval(rho_zz_1(k,1:nCells)) /= minval(rho_zz_1(k,1:nCells))) then
               call mpas_log_write('dcmip2012_1_2 needs rho_zz uniform in the horizontal (init case 14).', &
                                   messageType=MPAS_LOG_CRIT)
            end if
         end do
      end if

      allocate(rho_col(nVertLevels), vmer(nVertLevels), s_int(nVertLevels+1))
      rho_col(:) = rho_zz_1(1:nVertLevels,1)

      ! Step mean of cos(pi t/tau), times the amplitude of rho v
      amp = -a * w0 * pii * rho0 / (kc * ztop) &
            * tau / (pii * dt) * (sin(pii * (t0 + dt) / tau) - sin(pii * t0 / tau))

      ! Layer means of cos(pi z/ztop); (ztop/pi) sin(pi z/ztop) at the interfaces (flat terrain: column 1)
      do k = 1, nVertLevels
         vmer(k) = ztop / (pii * (zgrid(k+1,1) - zgrid(k,1))) &
                   * (sin(pii * zgrid(k+1,1) / ztop) - sin(pii * zgrid(k,1) / ztop))
      end do
      do k = 1, nVertLevels+1
         s_int(k) = ztop / pii * sin(pii * zgrid(k,1) / ztop)
      end do

      ! Streamfunction of the zonal solid-body flow at the vertices: u = -(1/a) dpsi/dlat
      allocate(psi(nVertices+1))
      psi(:) = 0.0_RKIND
      do iVtx = 1, nVertices
         psi(iVtx) = -u0 * a * sin(latVertex(iVtx))
      end do

      !
      ! Horizontal mass fluxes
      !
      allocate(m_edge(nEdges+1))
      m_edge(:) = 0.0_RKIND
      n_fallback = 0
      err_u = 0.0_RKIND
      max_u = 0.0_RKIND
      do iEdge = 1, nEdges
         v1 = verticesOnEdge(1,iEdge)
         v2 = verticesOnEdge(2,iEdge)
         uz_mid = u0 * cos(latEdge(iEdge)) * cos(angleEdge(iEdge))
         mer_mid = cos(latEdge(iEdge)) * sin(kc * latEdge(iEdge)) * sin(angleEdge(iEdge))
         if (v1 >= 1 .and. v1 <= nVertices .and. v2 >= 1 .and. v2 <= nVertices) then
            uz = -(psi(v2) - psi(v1)) / dvEdge(iEdge)
            p1 = (/ xVertex(v1), yVertex(v1), zVertex(v1) /)
            p2 = (/ xVertex(v2), yVertex(v2), zVertex(v2) /)
            p1 = p1 / norm2(p1)
            p2 = p2 / norm2(p2)
            nrm = (/ p2(2)*p1(3) - p2(3)*p1(2), p2(3)*p1(1) - p2(1)*p1(3), p2(1)*p1(2) - p2(2)*p1(1) /)
            nrm = nrm / norm2(nrm)                       ! n = (v2 x v1)/|v2 x v1|, constant along the arc
            arc = acos(max(-1.0_RKIND, min(1.0_RKIND, dot_product(p1, p2))))
            mer = 0.0_RKIND
            do is = 1, 3
               s = 0.5_RKIND * (1.0_RKIND + gl_x(is))
               pq = (sin((1.0_RKIND - s) * arc) * p1 + sin(s * arc) * p2) / sin(arc)
               lat = asin(max(-1.0_RKIND, min(1.0_RKIND, pq(3))))
               lon = atan2(pq(2), pq(1))
               mer = mer + gl_w(is) * cos(lat) * sin(kc * lat) &
                         * (-sin(lat)*cos(lon) * nrm(1) - sin(lat)*sin(lon) * nrm(2) + cos(lat) * nrm(3))
            end do
            err_u = max(err_u, abs(uz - uz_mid), abs(amp * (mer - mer_mid)))
            max_u = max(max_u, abs(uz_mid), abs(amp * mer_mid))
         else
            uz = uz_mid                              ! outer-halo edge; never updates an owned cell
            mer = mer_mid
            n_fallback = n_fallback + 1
         end if
         m_edge(iEdge) = mer
         do k = 1, nVertLevels
            ruAvg(k,iEdge) = rho_col(k) * uz + amp * vmer(k) * mer
            u_2(k,iEdge) = ruAvg(k,iEdge) / rho_col(k)          ! diagnostic only: step-mean wind
         end do
      end do
      ruAvg(:,nEdges+1) = 0.0_RKIND

      !
      ! Vertical mass flux from the discrete horizontal divergence of g(lat) j
      !
      wwAvg(:,:) = 0.0_RKIND
      err_w = 0.0_RKIND
      max_w = 0.0_RKIND
      do iCell = 1, nCells
         cdiv_exact = (-2.0_RKIND * sin(kc * latCell(iCell)) * sin(latCell(iCell)) &
                       + kc * cos(latCell(iCell)) * cos(kc * latCell(iCell))) / a
         if (all(edgesOnCell(1:nEdgesOnCell(iCell),iCell) <= nEdges)) then
            cdiv = 0.0_RKIND
            do j = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(j,iCell)
               cdiv = cdiv + edgesOnCell_sign(j,iCell) * m_edge(iEdge) * dvEdge(iEdge)
            end do
            cdiv = cdiv * invAreaCell(iCell)
            ! peak of rho w over z is at sin(pi z/ztop) = 1: |amp| ztop/pi |C| (kg m-2 s-1)
            err_w = max(err_w, abs(amp * ztop / pii * (cdiv - cdiv_exact)))
            max_w = max(max_w, abs(amp * ztop / pii * cdiv_exact))
         else
            cdiv = cdiv_exact                         ! outer-halo cell; never updated as owned
         end if
         wwAvg(2:nVertLevels,iCell) = -amp * cdiv * s_int(2:nVertLevels)
      end do

      rho_zz_2(:,:) = rho_zz_1(:,:)

      if (check) then
         ! 3D mass divergence of the owned cells, relative to its horizontal part. Owned cells only:
         ! a halo cell next to a fallback edge would fail this check without affecting the solution.
         res = 0.0_RKIND
         res_scale = 0.0_RKIND
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               hdiv = 0.0_RKIND
               do j = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(j,iCell)
                  hdiv = hdiv + edgesOnCell_sign(j,iCell) * ruAvg(k,iEdge) * dvEdge(iEdge) * invAreaCell(iCell)
               end do
               res = max(res, abs(hdiv + rdzw(k) * (wwAvg(k+1,iCell) - wwAvg(k,iCell))))
               res_scale = max(res_scale, abs(hdiv))
            end do
         end do
         call mpas_log_write('Kinematic fluxes (dcmip2012_1_2): max |u.n| = $r m/s, max |flux - point wind| = $r m/s,' &
                             //' max |rho w| = $r kg/m2/s, max |ww - cell-centre rho w| = $r, fallback edges = $i', &
                             realArgs=(/max_u, err_u, max_w, err_w/), intArgs=(/n_fallback/))
         call mpas_log_write('Kinematic fluxes (dcmip2012_1_2): max |3D mass divergence| = $r, relative $r', &
                             realArgs=(/res, res / max(res_scale, tiny(1.0_RKIND))/))
         if (err_u > 0.1_RKIND * max_u .or. err_w > 0.1_RKIND * max_w) then
            call mpas_log_write('Kinematic fluxes disagree with the analytic wind; check the edge orientation.', &
                                messageType=MPAS_LOG_CRIT)
         end if
         if (res > 1.0e-12_RKIND * res_scale) then
            call mpas_log_write('dcmip2012_1_2 mass flux is not divergence-free.', messageType=MPAS_LOG_CRIT)
         end if
      end if

      deallocate(psi, m_edge, vmer, s_int, rho_col)

   end subroutine kinematic_dcmip2012_1_2

end module atm_kinematic
