! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vertical_remap
!
!> \brief MPAS ocean vertical Lagrangian remapping
!> \author Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis
!> \date   July 2021
!> \details
!>  This module contains the vertical remapping routine.
!
!-----------------------------------------------------------------------

module ocn_regrid_isopycnal

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics
   use ocn_diagnostics_variables

   implicit none
   private
   save

   contains

!***********************************************************************
!
!  routine ocn_vertical_target_isopycnal
!
!> \brief   
!> \author  Carolyn Begeman
!> \date    2021
!> \details
!>  This routine 
!
!-----------------------------------------------------------------------

   subroutine ocn_vertical_target_isopycnal(ssh, mask, layerThicknessLag, &
                                            layerThicknessTarget, err)

      integer, dimension(:), intent(in) :: &
         mask ! only solve for target layers in masked area
      real (kind=RKIND), dimension(:), intent(in) :: &
         ssh ! ssh at current timestep
      real (kind=RKIND), dimension(:, :), intent(in) :: &
         layerThicknessLag ! layerThickness after the lagrangian step
      real (kind=RKIND), dimension(:, :), intent(inout) :: &
         layerThicknessTarget ! target layerThickness for remapping
      integer, intent(out) :: err !< Output: Error flag

      integer :: i, iCell, k, kmin, kmax, kLayer, nContours, nFit
      real (kind=RKIND) :: minLayerThickness = 10.0, maxLayerThickness = 100.0
      real (kind=RKIND) :: rhominAll = 1000.0, rhomaxAll = 1040.0, rhoinc = 0.1
      real (kind=RKIND) :: zinc = 1.0, m
      real (kind=RKIND), dimension(:), allocatable :: &
         rhoContours, rhoContourLoc, rhoFit, zFit, zTop, zMid

      ! potentialDensity imported from ocn_diagnostics_variables
      !call ocn_equation_of_state_density(statePool, meshPool, activeTracers, tracersSurfaceValue, &
      !                                   nCells, 1, 'absolute', potentialDensity, &
      !                                   err, timeLevelIn=timeLevel)
      nContours = (rhomaxAll-rhominAll)/rhoinc
      allocate(rhoContours(nContours))
      allocate(rhoContourLoc(nVertLevels))
      allocate(zTop(nVertLevels), zMid(nVertLevels))
      rhoContours = [(rhominAll + (i-1) * rhoinc, i=1,nContours)]
      do iCell = 1, nCellsAll
         if (mask(iCell)) then

            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)

            ! Start with piecewise linear construction
            nFit = FLOOR((bottomDepth(iCell)+ssh(iCell))/zinc)
            allocate(zFit(nFit),rhoFit(nFit))
            zFit(1) = -ssh(iCell)
            do k = 2,nFit
               zFit(k) = zFit(k-1) + zinc
            enddo
            zTop(kmin) = -ssh(iCell)
            do k = kmin+1,kmax+1
               zMid(k-1) = zTop(k-1) + layerThicknessLag(k-1,iCell)/2
               zTop(k) = zTop(k-1) + layerThicknessLag(k-1,iCell)
            enddo
            do k = 1, nFit
               kLayer = MINLOC(ABS(zMid(k)-zFit(k)))
               if (kLayer == kmax) then
                  rhoFit(k) = potentialDensity(k,iCell)
               else
                  m = (potentialDensity(klayer+1,iCell) - potentialDensity(klayer,iCell))/&
                      (zTop(kLayer+1) - zTop(kLayer))
                  rhoFit(k) = m*zFit(k) + (potentialDensity(kLayer,iCell) - m*zTop(kLayer))
               endif
            enddo
            
            ! Adjust min and max limits
            ! TODO work out syntax
            !do while (SUM(refLayerThickness(1:kmin-1,iCell)) > ssh(iCell))
            !   kmin = kmin+1 
            !enddo
            !do while (SUM(refLayerThickness(1:kmin,iCell)) <= &
            !          ssh(iCell)+maxLayerThickness)
            !   kmin = kmin-1 
            !enddo

            rhomin = MIN(potentialDensity(kmin:kmax,iCell))
            rhomax = MAX(potentialDensity(kmin:kmax,iCell))
            rhoContourLoc(kmin-1) = rhoContours(MINLOC(ABS(rhoContours-rhomin)))
            zTop(kmin) = -ssh(iCell)
            do k = kmin,kmax
               rhoContourLoc(k) = rhoContourLoc(k-1) + rhoinc
               zTop(k+1) = MAX(zTop(k)+minLayerThickness, &
                               zFit(MINLOC(rhoFit-rhoContourLoc)))
               zTop(k+1) = MIN(zTop(k)+maxLayerThickness, &
                               zTop(k+1))
               zTop(k+1) = MIN(bottomDepth(iCell), &
                               zTop(k+1))
               layerThicknessTarget(k,iCell) = zTop(k+1) - zTop(k)
               if (layerThicknessTarget(k,iCell) < minLayerThickness) then
                  kmax = k
                  exit
               endif
            enddo
            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)
            deallocate(zFit,rhoFit)
         endif
      enddo

   end subroutine
end module ocn_regrid_isopycnal
! vim: foldmethod=marker
