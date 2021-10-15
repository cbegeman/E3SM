! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_sediment
!
!> \brief MPAS ocean sediment
!> \author Zhendong Cao, Mathew Maltrud
!> \date   02/04/2021
!> \details
!>  This module contains routines for computing tracer forcing due to sediment
!
!-----------------------------------------------------------------------

module ocn_tracer_sediment

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_config
   use ocn_diagnostics
   use ocn_diagnostics_variables

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_sediment_compute, &
             ocn_tracer_sediment_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   real(kind=RKIND) :: rho_s, bottomDrag
   real(kind=RKIND), parameter :: rho_w = 1020, gravity = 9.81
   real(kind=RKIND), dimension(:), allocatable :: &
      wSettling, sedimentThickness

contains

!***********************************************************************
!
!  routine ocn_tracer_sediment_compute
!
!> \brief   computes a tracer tendency due to sediment
!> \author  Zhendong Cao, Mathew Maltrud
!> \date    02/04/2021
!> \details
!>  This routine computes a tracer tendency due to sediment
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sediment_compute(sedimentTracers, &
      nTracersSediment, dt, nCellsSolve, minLevelCell, &
      maxLevelCell, nVertLevels, layerThickness, w, normalVelocity, &
      sedimentTracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: dt

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         minLevelCell, maxLevelCell

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, normalVelocity

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         sedimentTracers

      ! scalars
      integer, intent(in) :: nTracersSediment, nCellsSolve, nVertLevels

      !
      ! two dimensional pointers
      !
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         w ! total diasurface vertical velocity

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         sedimentTracersTend

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! source/sink wants cm instead of m

      integer :: iCell, k, iTrac

      real (kind=RKIND), parameter :: &
         porosity = 0.5, &
         tau_crit = 0.1, &
         min_erosion_rate = 5e-4

      real (kind=RKIND) :: tau_mean
      
      real (kind=RKIND), dimension(:), allocatable :: &
         erosionFlux, settlingFlux, bedFlux, &
         sedimentThicknessTend
      
      write(*,*) 'Enter sediment_tracer_compute'

      if (.not. config_sediment_suspended_enabled) return

      err = 0

      call mpas_timer_start("sediment source-sink")

      iTrac = 1

      ! Vertical velocity of sediment is updated to include both vertical
      ! advection applied to all tracers and settling
      write(*,*) '-- add wSettling'
      do iCell=1,nCellsSolve
         ! w is the velocity through the top interface and wSettling applies to
         ! bottom interface
         do k=minLevelCell(iCell)+1, maxLevelCell(iCell)
            w(k,iCell) = w(k,iCell) - wSettling(iTrac)
         enddo  !  iLevel

      enddo  !  iCell

      call mpas_timer_stop("sediment source-sink")

      if (config_sediment_bedload_enabled) then

         call mpas_timer_start("sediment surface flux")
         allocate(sedimentThicknessTend(nCellsSolve))
         allocate(erosionFlux(nCellsSolve))
         allocate(settlingFlux(nCellsSolve))
         allocate(bedFlux(nCellsSolve))
         call ocn_diagnostic_solve_kineticEnergy(normalVelocity, kineticEnergyCell)

         !DWJ 08/05/2016: This loop needs OpenMP added to it
         ! Only the tendency at the bottom of the domain is computed
         write(*,*) '-- compute settlingFlux'
         do iCell = 1, nCellsSolve
            k = maxLevelCell(iCell)
            ! settling flux
            ! - maximum settling flux is that which settles all sediment in the bottom cell
            settlingFlux(iCell) = MIN(wSettling(iTrac) * &
                                      sedimentTracers(iTrac,k,iCell), &
                                      sedimentTracers(iTrac,k,iCell) * &
                                      layerThickness(k,iCell)/dt )
         enddo  !  iCell
            
         write(*,*) '-- compute erosionFlux'
         do iCell = 1, nCellsSolve
            k = maxLevelCell(iCell)
            ! surface erosion mass flux ero_flux in unit (kg m-2 s-1)
            ! - maximum erosion flux uses all sediment available
            
            ! compute shear stress [m^2 s^{-2}]
            ! TODO Make sure that tidal velocity is accounted for in kineticEnergyCell
            tau_mean = bottomDrag * kineticEnergyCell(k,iCell)
            erosionFlux(iCell) = MIN((1.0_RKIND - porosity) * min_erosion_rate * &
                                      MAX(tau_mean/tau_crit - 1.0_RKIND, 0.0_RKIND), &
                                      sedimentThickness(iCell)*rho_s/dt)
         enddo  !  iCell
         write(*,*) '-- compute bedFlux'
         do iCell = 1, nCellsSolve
            k = maxLevelCell(iCell)
            ! total flux to the bed
            bedFlux(iCell) = erosionFlux(iCell) - settlingFlux(iCell)
         enddo  !  iCell
         write(*,*) '-- compute sedimentTracersTend'
         do iCell = 1, nCellsSolve
            k = maxLevelCell(iCell)
            sedimentTracersTend(iTrac,k,iCell) = sedimentTracersTend(iTrac,k,iCell) + &
                                                 dt * bedFlux(iCell) / &
                                                 layerThickness(k,iCell)
         enddo  !  iCell
         write(*,*) 'bedFlux',bedFlux(1), 'settlingFlux', settlingFlux(1), 'erosionFlux', erosionFlux(1)
         write(*,*) 'sedimentTracersTend',sedimentTracersTend(1,1,1)

         !write(*,*) '-- compute sedimentThicknessTend'
         !do iCell = 1, nCellsSolve
         !   k = maxLevelCell(iCell)
         !   ! compute diagnostic related to bed fluxes
         !   ! TODO store sedimentThicknessTend as sedimentTracersBedTend
         !   sedimentThicknessTend(iCell) = dt * bedFlux(iCell)/ &
         !                                  (rho_s*(1.0_RKIND - porosity))
         !enddo  !  iCell
         !write(*,*) '-- compute sedimentThickness'
         !do iCell = 1, nCellsSolve
         !   k = maxLevelCell(iCell)
         !   sedimentThickness(iCell) = MAX(0.0_RKIND, &
         !                                sedimentThickness(iCell) &
         !                              + sedimentThicknessTend(iCell))

         !enddo  !  iCell
         deallocate(sedimentThicknessTend, bedFlux, erosionFlux, settlingFlux)
         call mpas_timer_stop("sediment surface flux")
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_sediment_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_sediment_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Zhendong Cao, Mathew Maltrud
!> \date    02/04/2021
!> \details
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sediment_init(domain, err)!{{{

!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        sedimentTracers

      ! scalars
      integer :: iTrac, nTracersSediment=1

      ! scalar pointers
      integer, pointer :: nVertLevels, nCells

      real (kind=RKIND), parameter :: &
         Vh = 1e-6_RKIND        ! kinematic viscosity
      real (kind=RKIND) :: &
         SD50, rho_R, ND50, dstar
      !
      ! get tracers pools
      !

      err = 0

      !
      ! Get tracer group so we can get the number of tracers in it
      !

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'sedimentTracers', sedimentTracers, 1)

      if (associated(sedimentTracers)) then

         sedimentTracers = config_sediment_suspended_initConcen
         SD50 = config_sediment_grain_diameter
         rho_s = config_sediment_grain_density

         ! sieve diameter (SD50) to nominal diameter (ND50), Raudkivi [1990]
         ND50 = SD50*1.1_RKIND
         ! sediment/water density comparison  -> specific graivity
         rho_R = rho_s/rho_w - 1! rhoS/rho0-1

         call mpas_pool_get_dimension(domain % blocklist % dimensions, &
                                      'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(domain % blocklist % dimensions, &
                                      'nCells', nCells)

         ! COMPUTE Sediment settling velocity [m/s] 
         allocate(wSettling(nTracersSediment))
         wSettling(:) = 0.0_RKIND
         
         iTrac = 1 ! for now, only one sediment class
         if (trim(config_sediment_settling_formulation) .ne. "None" .and. &
             trim(config_sediment_settling_formulation) .ne. "Goldstein-Coco" .and. &
             trim(config_sediment_settling_formulation) .ne. "Cheng" .and. &
             trim(config_sediment_settling_formulation) .ne. "Soulsby") then
            call mpas_log_write( &
               "Unknown value for config_sediment_settling_formulation, " // &
               "supported values are:" // &
               "         None, Goldstein-Coco, Cheng, or Soulsby", &
               MPAS_LOG_CRIT)
            err = 1
            return
         endif
         if (config_sediment_settling_formulation == 'Goldstein-Coco') then
            ! Goldstein & Coco [2013]
            wSettling(iTrac) = (37.8_RKIND*ND50*rho_R*(1 + 100.0_RKIND*ND50))/ &
                                 (0.383_RKIND +                                 & 
                                  1E4_RKIND*rho_R*Vh +                          &
                                  1E2_RKIND*ND50*rho_R**2)
         elseif (config_sediment_settling_formulation == 'Soulsby') then
            wSettling(iTrac) = 10.36_RKIND*Vh/SD50 *                           &
                                (sqrt(1.0_RKIND +                               &
                                      0.156_RKIND*rho_R*gravity*SD50**3.0_RKIND/&
                                      (16.0_RKIND*Vh**2.0) ) - 1)
         elseif (config_sediment_settling_formulation == 'Cheng') then
            dstar = SD50 * (rho_R*gravity/Vh**2.0_RKIND)**0.333_RKIND
            wSettling(iTrac) = Vh/SD50 * (sqrt(25.0_RKIND+1.2_RKIND*dstar**2) -&
                                           5.0_RKIND)**1.5_RKIND
         endif

         allocate(sedimentThickness(nCells))
         sedimentThickness(:) = config_sediment_initial_bed_thickness

      end if  !  associated(sedimentTracers)

      if (config_use_implicit_bottom_drag) then
         bottomDrag = config_implicit_bottom_drag_coeff
      else
         bottomDrag = config_explicit_bottom_drag_coeff
      endif

   end subroutine ocn_tracer_sediment_init!}}}

!***********************************************************************

end module ocn_tracer_sediment

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
