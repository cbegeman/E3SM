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

! ZHENDONG: need turn those on or define those here later ...   
   use SEDIMENT_mod
   use SEDIMENT_parms
   !use BGC_mod
   !use BGC_parms

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
             ocn_tracer_sediment_surface_flux_compute,  &
             ocn_tracer_sediment_init

   integer, public:: &
      nCells

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  name the necessary sediment derived types
!  all of these are defined in sediment_mod
!-----------------------------------------------------------------------

  type(SEDIMENT_indices_type)    , public :: SEDIMENT_indices
  type(SEDIMENT_input_type)      , public :: SEDIMENT_input
  type(SEDIMENT_output_type)     , public :: SEDIMENT_output
  type(SEDIMENT_diagnostics_type), public :: SEDIMENT_diagnostic_fields

! hold indices in tracer pool corresponding to each tracer array
  type(SEDIMENT_indices_type), public :: sedimentIndices

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_sediment_compute
!
!> \brief   computes a tracer tendency due to sediment
!> \author  Zhendong Cao, Mathew Maltrud
!> \date    02/04/2021
!> \details
!>  This routine only places bounds on sediment concentration
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sediment_compute(sedimentTracers, nTracersSediment, &
      nCellsSolve, minLevelCell, maxLevelCell, nVertLevels, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         minLevelCell, maxLevelCell

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

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         sedimentTracers

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

      integer :: iCell, iLevel, iTracer

      call mpas_timer_start("sediment source-sink")

      err = 0

      do iCell=1,nCellsSolve
         SEDIMENT_input%number_of_active_levels(iCell) = maxLevelCell(iCell) - minLevelCell(iCell) + 1
         do iTracer=1,nTracersSediment
         SEDIMENT_input%SEDIMENT_tracers(:,iCell,iTracer) = &
            sedimentTracers(iTracer,:,iCell)
         enddo
         call SEDIMENT_SourceSink(SEDIMENT_indices, SEDIMENT_input, &
                                  SEDIMENT_output, SEDIMENT_diagnostic_fields,           &
                                  nVertLevels, nCellsSolve)
         do iTracer=1,nTracersSediment
         sedimentTracers(iTracer,:,iCell) = &
            SEDIMENT_output%SEDIMENT_tracers(:,iCell,iTracer)
         enddo

      enddo  !  iCell

      call mpas_timer_stop("sediment source-sink")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_sediment_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_sediment_surface_flux_compute
!
!> \brief   computes a tracer tendency due to sediment
!> \author  Zhendong Cao, Mathew Maltrud
!> \date    02/04/2021
!> \details
!>  This routine computes a tracer tendency due to sediment
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sediment_surface_flux_compute(sedimentTracers, &
                                       nCellsSolve, zMid, sedimentSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         sedimentSurfaceFlux

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         sedimentTracers

      ! scalars
      integer, intent(in) :: nCellsSolve


      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      integer :: iCell

      call mpas_timer_start("sediment surface flux")

      err = 0

      ! fluxes are zero

      !DWJ 08/05/2016: This loop needs OpenMP added to it
      do iCell = 1, nCellsSolve

      enddo  !  iCell

      call mpas_timer_stop("sediment surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_sediment_surface_flux_compute!}}}

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

   subroutine ocn_tracer_sediment_init(domain,err)!{{{

!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        sedimentTracers

      ! scalars
      integer :: nCells

      ! scalar pointers
      integer, pointer :: nVertLevels, index_dummy

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

      ! make sure sediment is turned on

      if (associated(sedimentTracers)) then

      !
      ! pull nVertLevels out of the mesh structure
      !

         call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nCells', nCells)

!-----------------------------------------------------------------------
!  initialize sediment parameters
!-----------------------------------------------------------------------

         allocate( SEDIMENT_indices%short_name(SEDIMENT_tracer_cnt) )
         allocate( SEDIMENT_indices%long_name(SEDIMENT_tracer_cnt) )
         allocate( SEDIMENT_indices%units(SEDIMENT_tracer_cnt) )


!-----------------------------------------------------------------------
!  sets most of sediment parameters
!  sets namelist defaults
!-----------------------------------------------------------------------

         call SEDIMENT_parms_init

         SEDIMENT_indices%sedConcen_ind      = 1

         call mpas_pool_get_dimension(tracersPool, 'index_sedConcen',  index_dummy)
           sedimentIndices%sedConcen_ind  = index_dummy
! sediment_init sets short and long names, units in SEDIMENT_indices
         call SEDIMENT_init(SEDIMENT_indices)

!NOTES:

!also check short_name with mpas variable name

!-----------------------------------------------------------------------
!  allocate input, forcing, diagnostic arrays
!-----------------------------------------------------------------------

         allocate ( SEDIMENT_input%SEDIMENT_tracers(nVertLevels, nCells, SEDIMENT_tracer_cnt) )
         allocate ( SEDIMENT_input%cell_thickness(nVertLevels, nCells) )
         allocate ( SEDIMENT_input%number_of_active_levels(nCells) )

         allocate ( SEDIMENT_output%SEDIMENT_tendencies(nVertLevels, nCells, SEDIMENT_tracer_cnt) )
         allocate ( SEDIMENT_output%SEDIMENT_tracers(nVertLevels, nCells, SEDIMENT_tracer_cnt) )

    !---------------------------------------------------------------------------
    !   allocate diagnostic output fields
    !---------------------------------------------------------------------------
         SEDIMENT_input%number_of_active_tracers = config_sediment_size_classes


      end if  !  associated(sedimentTracers)



   !--------------------------------------------------------------------

   end subroutine ocn_tracer_sediment_init!}}}

!***********************************************************************

end module ocn_tracer_sediment

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
