! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vertical_regrid
!
!> \brief MPAS ocean vertical regridding
!> \author Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis
!> \date   July 2021
!> \details
!>  This module contains the vertical regridding routine, used for
!>  vertical Lagrangian remapping.
!
!-----------------------------------------------------------------------

module ocn_vertical_regrid

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_log

   use ocn_constants
   use ocn_config
   use ocn_mesh

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   character (len=100) :: log_string
   real (kind=RKIND) :: minLayerThickness, maxLayerThicknessFactor
   real (kind=RKIND) :: tau, mu ! config options for regridding
   logical, parameter :: zstar_coordinate = .true.
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vert_regrid
   public :: ocn_vert_regrid_init

   contains

!***********************************************************************
!
!  routine ocn_vert_regrid
!
!> \brief   Regridding for vertical Lagrangian remapping
!> \author  Carolyn Begeman
!> \date    2021
!> \details
!>  This routine determines the layerThickness to remap to subject to
!>  constraints
!
!-----------------------------------------------------------------------

   subroutine ocn_vert_regrid(mask, restingThickness, layerThicknessLag, &
                              layerThicknessTarget, err)

      integer, dimension(:), intent(in) :: &
         mask ! only adjust target layers in masked area (mask=1)
      real (kind=RKIND), dimension(:, :), intent(in) :: &
         layerThicknessLag ! layerThickness after the lagrangian step
      real (kind=RKIND), dimension(:, :), intent(inout) :: &
         restingThickness  ! updated if adding new layer
      real (kind=RKIND), dimension(:, :), intent(out) :: &
         layerThicknessTarget ! adjusted target layerThickness for remapping
      integer, intent(out) :: err !< Output: Error flag

      integer :: iCell, k, kmin, kmax

      real (kind=RKIND), dimension(:), allocatable :: &
            layerScalingFactor, totalRestingThickness, totalThickness
      real (kind=RKIND), dimension(:, :), allocatable :: &
            zLag, zTarget

      allocate(totalThickness(nCellsAll), &
               totalRestingThickness(nCellsAll), &
               layerScalingFactor(nCellsAll))

      ! --- Step 1: calculate z-star layer locations ---
      layerThicknessTarget = 0.0_RKIND

      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1,nCellsAll
         totalThickness(iCell) = 0.0_RKIND
         totalRestingThickness(iCell) = 0.0_RKIND
         do k = maxLevelCell(iCell), minLevelCell(iCell), -1
            totalRestingThickness(iCell) = totalRestingThickness(iCell) &
                                    + restingThickness(k,iCell)
            totalThickness(iCell) = totalThickness(iCell) &
                                    + layerThicknessLag(k,iCell)
         end do
         layerScalingFactor(iCell) = totalThickness(iCell)/ &
                                     totalRestingThickness(iCell)

         !if (iCell == 7) then
         !   write(log_string,*) 'Regrid S1: H,Hrest=',&
         !      totalThickness(iCell),totalRestingThickness(iCell)
         !   call mpas_log_write(log_string)
         !   do k = maxLevelCell(iCell), minLevelCell(iCell), -1
         !   write(log_string,*) 'Regrid S1: k,hrest=',&
         !      k,restingThickness(k,iCell)
         !   call mpas_log_write(log_string)
         !   end do
         !endif
         do k = maxLevelCell(iCell), minLevelCell(iCell), -1
            layerThicknessTarget(k, iCell) = restingThickness(k, iCell) * &
                                             layerScalingFactor(iCell)
            !if (iCell == 7) then
            !   write(log_string,*) 'Regrid S1: k,h(k)=',k,layerThicknessTarget(k,iCell)
            !   call mpas_log_write(log_string)
            !endif
         end do

      end do
      !$omp end do
      !$omp end parallel

      ! --- Step 2: Add/drop layers ---

      if (config_use_min_max_thickness) &
         call ocn_vert_regrid_minmax_thickness(mask, restingThickness, &
            layerThicknessTarget, err)

   end subroutine ocn_vert_regrid

!***********************************************************************
!
!  routine ocn_vert_regrid_minmax_thickness
!
!> \brief   
!> \author  Carolyn Begeman
!> \date    2021
!> \details
!>  This routine determines the layerThickness to remap to subject to
!>  constraints
!
!-----------------------------------------------------------------------

   subroutine ocn_vert_regrid_minmax_thickness( &
         mask, restingThickness, layerThicknessTarget, err)

      integer, dimension(:), intent(in) :: &
         mask ! only adjust target layers in masked area (mask=1)
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         restingThickness, &  ! If needed, fill in restingThickness for newly
                              ! created layers
         layerThicknessTarget ! target layerThickness for remapping
      integer, intent(out) :: err !< Output: Error flag

      logical :: adjust_kmin, adjust_kmax
      integer, parameter :: minLayers = 3
      integer :: iCell, iEdge, iNeighbor, k, kmin, kmax, nEdges
      real (kind=RKIND), dimension(:), allocatable :: &
            layerScalingFactor, totalRestingThickness, totalThickness
      real (kind=RKIND), dimension(:,:), allocatable :: &
         layerThicknessFactor

      err = 0

      allocate(layerThicknessFactor(nVertLevels,1))
      allocate(totalThickness(nCellsAll), &
               totalRestingThickness(nCellsAll), &
               layerScalingFactor(nCellsAll))
      
      !$omp parallel
      !$omp do schedule(runtime) private(k,kmin,kmax,adjust_kmin,adjust_kmax,layerThicknessFactor)
      do iCell = 1, nCellsAll
         if (mask(iCell)==0) cycle
         adjust_kmin = .false.
         adjust_kmax = .false.

         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)

         ! Evaluate layerThicknessFactor
         layerThicknessFactor(:,1) = 0.0_RKIND
         do k = kmin, kmax
            layerThicknessFactor(k,1) = layerThicknessTarget(k,iCell)/ &
                                        restingThickness(k,iCell)
         enddo
 
         ! Add a layer
         if (MAXVAL(layerThicknessFactor(kmin:kmax,1)) > &
             maxLayerThicknessFactor .and. zstar_coordinate) then
            ! Add layer from the bottom only if cannot add from top
            if (kmin == 1 .and. kmax < nVertLevels) then
               adjust_kmax = .true.
               kmax = kmax + 1
               write(log_string,*) 'add bottom layer iCell,mlc',iCell,maxLevelCell(iCell)
               call mpas_log_write(log_string)
               if (restingThickness(kmax,iCell) <= 0.0_RKIND) &
                  restingThickness(kmax,iCell) = restingThickness(kmax-1,iCell)
               !nEdges = nEdgesOnCell(iCell)
               !do iEdge = 1,nEdges
               !   iNeighbor = cellsOnCell(iEdge,iCell)
               !   maxLevelEdgeBot(edgesOnCell(iEdge,iCell)) = &
               !      MAX(kmax, maxLevelCell(iNeighbor))
               !   maxLevelEdgeTop(edgesOnCell(iEdge,iCell)) = &
               !      MIN(kmax, maxLevelCell(iNeighbor))
               !enddo

            ! Add a layer from the top
            elseif (minLevelCell(iCell) > 1) then
               adjust_kmin = .true.
               kmin = kmin - 1
               write(log_string,*) 'add top layer iCell,mlc',iCell,minLevelCell(iCell)
               call mpas_log_write(log_string)
               !nEdges = nEdgesOnCell(iCell)
               !do iEdge = 1,nEdges
               !   iNeighbor = cellsOnCell(iEdge,iCell)
               !   minLevelEdgeBot(edgesOnCell(iEdge,iCell)) = &
               !      MAX(kmin, minLevelCell(iNeighbor))
               !   minLevelEdgeTop(edgesOnCell(iEdge,iCell)) = &
               !      MIN(kmin, minLevelCell(iNeighbor))
               !enddo
               if (restingThickness(kmin,iCell) <= 0.0_RKIND) &
                  restingThickness(kmin,iCell) = restingThickness(kmin+1,iCell)

            endif
            exit
         endif 
         if (MINVAL(layerThicknessTarget(kmin:kmax,iCell)) < minLayerThickness &
             .and. kmax - kmin + 1 > minLayers &
             .and. iCell == 4) then

            adjust_kmax=.true.
            ! Don't drop a layer from the top for ice shelf cavity
            ! to preserve terrain-following coordinate

            ! Drop a layer from the bottom
            write(log_string,*) 'drop bottom layer iCell,kmax',iCell,kmax
            call mpas_log_write(log_string)
            kmax = kmax - 1
            !nEdges = nEdgesOnCell(iCell)
            !do iEdge = 1,nEdges
            !   iNeighbor = cellsOnCell(iEdge,iCell)
            !      maxLevelEdgeBot(edgesOnCell(iEdge,iCell)) = &
            !         MAX(kmax, maxLevelCell(iNeighbor))
            !      maxLevelEdgeTop(edgesOnCell(iEdge,iCell)) = &
            !         MIN(kmax, maxLevelCell(iNeighbor))
            !enddo
            !exit
         endif 

         if (adjust_kmin .or. adjust_kmax) then
            ! -- Recompute zstar locations --

            ! Note: layerThicknessTarget has not been updated yet
            ! so use the original min,maxLevelCell limits
            totalThickness(iCell) = 0.0_RKIND
            do k = maxLevelCell(iCell), minLevelCell(iCell), -1
               totalThickness(iCell) = totalThickness(iCell) &
                                       + layerThicknessTarget(k,iCell)
            end do

            ! This totalRestingThickness may have changed with kmin,kmax
            totalRestingThickness(iCell) = 0.0_RKIND
            do k = kmax, kmin, -1
               totalRestingThickness(iCell) = totalRestingThickness(iCell) &
                                       + restingThickness(k,iCell)
            end do

            layerScalingFactor(iCell) = totalThickness(iCell)/ &
                                        totalRestingThickness(iCell)
            !if (iCell == 7) then
            !   write(log_string,*) 'Minmax end: H,Hrest=',&
            !      totalThickness(iCell),totalRestingThickness(iCell)
            !   call mpas_log_write(log_string)
            !   write(log_string,*) 'Minmax end: hrest(kmax),hrest(kmin)=',&
            !      restingThickness(maxLevelCell(iCell),iCell),restingThickness(minLevelCell(iCell),iCell)
            !   call mpas_log_write(log_string)
            !endif

            layerThicknessTarget(:,iCell) = 0.0_RKIND
            do k = maxLevelCell(iCell), minLevelCell(iCell), -1
               layerThicknessTarget(k, iCell) = restingThickness(k, iCell) * &
                                                layerScalingFactor(iCell)
               !if (iCell == 7) then
               !   write(log_string,*) 'Minmax end: k,h(k)=',k,layerThicknessTarget(k,iCell)
               !   call mpas_log_write(log_string)
               !endif
            end do

            minLevelCell(iCell) = kmin
            maxLevelCell(iCell) = kmax
         endif
      enddo
      !$omp end do
      !$omp end parallel

   end subroutine ocn_vert_regrid_minmax_thickness

!***********************************************************************
!
!  routine ocn_vert_regrid_init
!
!> \brief   Initializes ocean vertical regridding
!> \author  Carolyn Begeman
!> \date    July 2021
!> \details
!>  This routine initializes parameters required for vertical Lagrangian
!>  regridding
!
!-----------------------------------------------------------------------

   subroutine ocn_vert_regrid_init(err)

      integer, intent(out) :: err !< Output: Error flag

      err = 0
 
      ! Options for regridding
      tau = config_vert_remap_time_relax_factor
      mu = config_vert_remap_horizontal_diffusivity
      ! consider checking that these options bound the restingThickness
      if (config_use_min_max_thickness) then
         minLayerThickness = config_min_thickness
         maxLayerThicknessFactor = config_max_thickness_factor
      endif

   end subroutine ocn_vert_regrid_init

end module ocn_vertical_regrid
! vim: foldmethod=marker
