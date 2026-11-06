! Copyright (c) 2026, University of Michigan
!
! Initial conditions for the kinematic transport tests (config_init_case = 14).
!
module init_atm_kinematic_tests

   use mpas_kind_types, only : RKIND, StrKIND
   use mpas_derived_types, only : mpas_pool_type, MPAS_LOG_CRIT
   use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_array, mpas_pool_get_dimension
   use mpas_constants, only : pii, gravity, rgas, cp, omega
   use mpas_log, only : mpas_log_write
   use atm_advection, only : atm_initialize_advection_rk, atm_initialize_deformation_weights

   implicit none

   private
   public :: init_atm_case_kinematic, init_atm_passive_tracers


contains


   !-----------------------------------------------------------------------
   !  routine init_atm_case_kinematic
   !
   !> \brief Initial state for the kinematic transport tests
   !> \details
   !>  Flat terrain; nVertLevels uniform layers up to config_ztop (height
   !>  coordinate, zz = 1); isothermal, hydrostatic, dry atmosphere at rest
   !>  with T0 = 300 K and p = p0 exp(-g z / (Rd T0)), as in DCMIP-2012
   !>  tests 1-1 to 1-3 (Kent et al. 2014). The base state equals the full
   !>  state, so all perturbation fields are zero. Winds are zero: in
   !>  kinematic mode the model prescribes the mass fluxes from the first
   !>  step. The passive tracers are set by config_kinematic_case.
   !>  The state is not in discrete hydrostatic balance; it is meant only
   !>  for runs with config_kinematic_case /= 'off'.
   !-----------------------------------------------------------------------
   subroutine init_atm_case_kinematic(mesh, nCells, nVertLevels, state, diag, configs)

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: configs

      real (kind=RKIND), parameter :: p0 = 1.0e+05_RKIND, t0 = 300.0_RKIND

      character (len=StrKIND), pointer :: config_kinematic_case, config_interface_projection
      integer, pointer :: config_num_passive_tracers
      real (kind=RKIND), pointer :: config_ztop

      integer, pointer :: nEdges, nVertices, maxEdges
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: sphere_radius, nominalMinDc, cf1, cf2, cf3

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge, latVertex, lonVertex, fEdge, fVertex
      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp, surface_pressure
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx, dss
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3, scalars
      real (kind=RKIND), dimension(:,:), pointer :: pb, rb, rr, tb, rtb, rt, p, ppb, pp, rho, theta
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, theta_m, u, w, ru, rw, v

      real (kind=RKIND), dimension(nVertLevels+1) :: zw
      real (kind=RKIND), dimension(nVertLevels) :: dzw
      real (kind=RKIND) :: dz, zmid, pres, exner_k, cof1, cof2
      integer :: iCell, iEdge, iVtx, k, nz1

      nz1 = nVertLevels

      call mpas_pool_get_config(configs, 'config_kinematic_case', config_kinematic_case)
      call mpas_pool_get_config(configs, 'config_interface_projection', config_interface_projection)
      call mpas_pool_get_config(configs, 'config_num_passive_tracers', config_num_passive_tracers)
      call mpas_pool_get_config(configs, 'config_ztop', config_ztop)

      if (nz1 < 3) then
         call mpas_log_write('Init case 14 needs config_nvertlevels >= 3.', messageType=MPAS_LOG_CRIT)
      end if
      if (config_num_passive_tracers < 1) then
         call mpas_log_write('Init case 14 needs config_num_passive_tracers >= 1.', messageType=MPAS_LOG_CRIT)
      end if

      !
      ! Scale distances and areas from the unit sphere to sphere_radius (as in init_atm_case_jw)
      !
      call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
      call mpas_pool_get_array(mesh, 'xCell', xCell)
      call mpas_pool_get_array(mesh, 'yCell', yCell)
      call mpas_pool_get_array(mesh, 'zCell', zCell)
      call mpas_pool_get_array(mesh, 'xEdge', xEdge)
      call mpas_pool_get_array(mesh, 'yEdge', yEdge)
      call mpas_pool_get_array(mesh, 'zEdge', zEdge)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'nominalMinDc', nominalMinDc)

      xCell(:) = xCell(:) * sphere_radius
      yCell(:) = yCell(:) * sphere_radius
      zCell(:) = zCell(:) * sphere_radius
      xVertex(:) = xVertex(:) * sphere_radius
      yVertex(:) = yVertex(:) * sphere_radius
      zVertex(:) = zVertex(:) * sphere_radius
      xEdge(:) = xEdge(:) * sphere_radius
      yEdge(:) = yEdge(:) * sphere_radius
      zEdge(:) = zEdge(:) * sphere_radius
      dvEdge(:) = dvEdge(:) * sphere_radius
      dcEdge(:) = dcEdge(:) * sphere_radius
      areaCell(:) = areaCell(:) * sphere_radius**2.0
      areaTriangle(:) = areaTriangle(:) * sphere_radius**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * sphere_radius**2.0
      nominalMinDc = nominalMinDc * sphere_radius

      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius)
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

      !
      ! Coriolis parameter (unused in kinematic mode, set for completeness)
      !
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call mpas_pool_get_array(mesh, 'latVertex', latVertex)
      call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)
      do iEdge = 1, nEdges
         fEdge(iEdge) = 2.0_RKIND * omega * sin(latEdge(iEdge))
      end do
      do iVtx = 1, nVertices
         fVertex(iVtx) = 2.0_RKIND * omega * sin(latVertex(iVtx))
      end do

      !
      ! Vertical grid: flat terrain, uniform layers, height coordinate
      !
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'hx', hx)
      call mpas_pool_get_array(mesh, 'dss', dss)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

      dz = config_ztop / real(nz1, RKIND)
      do k = 1, nz1+1
         zw(k) = real(k-1, RKIND) * dz
      end do
      do k = 1, nz1
         dzw(k) = zw(k+1) - zw(k)
         rdzw(k) = 1.0_RKIND / dzw(k)
      end do
      dzu(1) = 0.0_RKIND
      rdzu(1) = 0.0_RKIND
      do k = 2, nz1
         dzu(k) = 0.5_RKIND * (dzw(k) + dzw(k-1))
         rdzu(k) = 1.0_RKIND / dzu(k)
      end do

      fzm(:) = 0.0_RKIND
      fzp(:) = 0.0_RKIND
      if (trim(config_interface_projection) == 'linear_interpolation') then
         do k = 2, nz1
            fzp(k) = 0.5_RKIND * dzw(k  ) / dzu(k)
            fzm(k) = 0.5_RKIND * dzw(k-1) / dzu(k)
         end do
      else if (trim(config_interface_projection) == 'layer_integral') then
         do k = 2, nz1
            fzm(k) = 0.5_RKIND * dzw(k  ) / dzu(k)
            fzp(k) = 0.5_RKIND * dzw(k-1) / dzu(k)
         end do
      else
         call mpas_log_write('only linear_interpolation or layer_integral are supported', messageType=MPAS_LOG_CRIT)
      end if

      cof1 = (2.0_RKIND*dzu(2) + dzu(3)) / (dzu(2) + dzu(3)) * dzw(1) / dzu(2)
      cof2 = dzu(2) / (dzu(2) + dzu(3)) * dzw(1) / dzu(3)
      cf1 = fzp(2) + cof1
      cf2 = fzm(2) - cof1 - cof2
      cf3 = cof2

      hx(:,:) = 0.0_RKIND
      zxu(:,:) = 0.0_RKIND
      dss(:,:) = 0.0_RKIND
      zb(:,:,:) = 0.0_RKIND        ! flat terrain: z at the edge equals z at the cells
      zb3(:,:,:) = 0.0_RKIND
      do iCell = 1, nCells
         zgrid(:,iCell) = zw(:)
         zz(:,iCell) = 1.0_RKIND
      end do

      !
      ! Isothermal hydrostatic state; base state = full state
      !
      call mpas_pool_get_array(diag, 'exner_base', pb)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'theta_base', tb)
      call mpas_pool_get_array(diag, 'rtheta_base', rtb)
      call mpas_pool_get_array(diag, 'rtheta_p', rt)
      call mpas_pool_get_array(diag, 'exner', p)
      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz)
      call mpas_pool_get_array(state, 'theta_m', theta_m)
      call mpas_pool_get_array(state, 'scalars', scalars)
      call mpas_pool_get_array(state, 'u', u)
      call mpas_pool_get_array(state, 'w', w)

      scalars(:,:,:) = 0.0_RKIND
      do iCell = 1, nCells
         do k = 1, nz1
            zmid = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
            pres = p0 * exp(-gravity * zmid / (rgas * t0))
            exner_k = (pres / p0)**(rgas / cp)
            ppb(k,iCell) = pres
            pp(k,iCell) = 0.0_RKIND
            pb(k,iCell) = exner_k
            p(k,iCell) = exner_k
            rb(k,iCell) = pres / (rgas * t0 * zz(k,iCell))
            rr(k,iCell) = 0.0_RKIND
            tb(k,iCell) = t0 / exner_k
            rtb(k,iCell) = rb(k,iCell) * tb(k,iCell)
            rt(k,iCell) = 0.0_RKIND
            rho_zz(k,iCell) = rb(k,iCell)
            theta_m(k,iCell) = tb(k,iCell)
         end do
      end do

      u(:,:) = 0.0_RKIND
      ru(:,:) = 0.0_RKIND
      v(:,:) = 0.0_RKIND
      w(:,:) = 0.0_RKIND
      rw(:,:) = 0.0_RKIND

      ! surface_pressure, rho and theta are set by decouple_variables, called after this routine

      !
      ! Tracers
      !
      select case (trim(config_kinematic_case))
      case ('lauritzen_nondiv')
         call init_atm_passive_tracers(mesh, state, configs, nCells)
      case default
         call mpas_log_write('Init case 14: unknown config_kinematic_case '''//trim(config_kinematic_case)//'''', &
                             messageType=MPAS_LOG_CRIT)
      end select

   end subroutine init_atm_case_kinematic


   !-----------------------------------------------------------------------
   !  routine init_atm_passive_tracers
   !
   !> \brief Initialize passive tracers with the standard transport-test shapes
   !> \details
   !>  Tracer n is set to shape mod(n-1,5), so that tracers n and n+5 are identical:
   !>    0: uniform value 1 (consistency with the mass field)
   !>    1: Gaussian hills      (Lauritzen et al. 2012, GMD, eqs. 7-9)
   !>    2: cosine bells        (eqs. 10-11)
   !>    3: slotted cylinders   (eq. 12)
   !>    4: correlated cosine bells, -0.8 * (cosine bells)**2 + 0.9 (eqs. 13-15);
   !>       with shape 2 it forms the tracer pair of the mixing diagnostics
   !>  The fields are uniform in the vertical. Each shape is a pair of
   !>  distributions centred at (5pi/6, 0) and (7pi/6, 0).
   !-----------------------------------------------------------------------
   subroutine init_atm_passive_tracers(mesh, state, configs, nCells)

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: nCells

      integer, pointer :: config_num_passive_tracers
      real (kind=RKIND), dimension(:,:,:), pointer :: passive_tracers
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell

      real (kind=RKIND), parameter :: radius = 0.5_RKIND   ! r/R for bells and cylinders
      real (kind=RKIND), dimension(2) :: lon_c, lat_c, x_c, y_c, z_c
      real (kind=RKIND), dimension(0:4) :: tracer_shape
      real (kind=RKIND) :: x, y, z, r_i, dlon
      integer :: iCell, itracer, i

      call mpas_pool_get_config(configs, 'config_num_passive_tracers', config_num_passive_tracers)
      if (config_num_passive_tracers <= 0) return

      call mpas_pool_get_array(state, 'passive_tracers', passive_tracers)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)

      lon_c(:) = (/ 5.0_RKIND*pii/6.0_RKIND, 7.0_RKIND*pii/6.0_RKIND /)
      lat_c(:) = 0.0_RKIND
      do i = 1, 2
         x_c(i) = cos(lat_c(i)) * cos(lon_c(i))
         y_c(i) = cos(lat_c(i)) * sin(lon_c(i))
         z_c(i) = sin(lat_c(i))
      end do

      do iCell = 1, nCells
         x = cos(latCell(iCell)) * cos(lonCell(iCell))
         y = cos(latCell(iCell)) * sin(lonCell(iCell))
         z = sin(latCell(iCell))

         tracer_shape(0) = 1.0_RKIND
         tracer_shape(1) = 0.0_RKIND
         tracer_shape(2) = 0.1_RKIND
         tracer_shape(3) = 0.1_RKIND

         do i = 1, 2
            ! Gaussian hills
            tracer_shape(1) = tracer_shape(1) &
                              + 0.95_RKIND * exp(-5.0_RKIND * ((x-x_c(i))**2 + (y-y_c(i))**2 + (z-z_c(i))**2))

            ! great-circle distance to centre i, and longitude offset wrapped to [-pi, pi)
            r_i = acos(max(-1.0_RKIND, min(1.0_RKIND, sin(lat_c(i))*sin(latCell(iCell)) &
                                         + cos(lat_c(i))*cos(latCell(iCell))*cos(lonCell(iCell)-lon_c(i)))))
            dlon = modulo(lonCell(iCell) - lon_c(i) + pii, 2.0_RKIND*pii) - pii

            ! cosine bells
            if (r_i < radius) then
               tracer_shape(2) = 0.1_RKIND + 0.9_RKIND * 0.5_RKIND * (1.0_RKIND + cos(pii * r_i / radius))
            end if

            ! slotted cylinders: slot open to the north on cylinder 1, to the south on cylinder 2
            if (r_i <= radius) then
               if (abs(dlon) >= radius/6.0_RKIND) then
                  tracer_shape(3) = 1.0_RKIND
               else if (i == 1 .and. latCell(iCell) - lat_c(1) < -5.0_RKIND/12.0_RKIND*radius) then
                  tracer_shape(3) = 1.0_RKIND
               else if (i == 2 .and. latCell(iCell) - lat_c(2) >  5.0_RKIND/12.0_RKIND*radius) then
                  tracer_shape(3) = 1.0_RKIND
               end if
            end if
         end do

         ! correlated cosine bells
         tracer_shape(4) = -0.8_RKIND * tracer_shape(2)**2 + 0.9_RKIND

         do itracer = 1, config_num_passive_tracers
            passive_tracers(itracer,:,iCell) = tracer_shape(mod(itracer-1, 5))
         end do
      end do

   end subroutine init_atm_passive_tracers

end module init_atm_kinematic_tests
