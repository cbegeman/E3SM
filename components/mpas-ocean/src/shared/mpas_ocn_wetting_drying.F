! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_wetting_drying
!
!> \brief MPAS ocean wetting and drying
!> \author Phillip J. Wolfram
!> \date   03/09/2018
!> \details
!>  This module contains the wetting and drying algorithms
!>  to prevent cell thicknesses from becoming too small.
!>
!
!-----------------------------------------------------------------------

module ocn_wetting_drying

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_io_units

   use ocn_constants
   use ocn_config
   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_gm
   use ocn_mesh
   use ocn_subgrid
   use ocn_vel_pressure_grad

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

   public :: ocn_wetting_drying_verify, ocn_prevent_drying_rk4, &
             ocn_wetting_drying_barotropic_factor

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   real(kind=RKIND), parameter :: eps = 1.0e-12_RKIND


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_wetting_drying_verify
!
!> \brief   Verifies that cells are not too dry.
!> \author  Phillip J. Wolfram
!> \date    03/09/2018
!> \details
!>  This routine checks that the minimum thickness in a cell is not
!>  too small.
!
!-----------------------------------------------------------------------

   subroutine ocn_wetting_drying_verify( block , minHeight, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (block_type), intent(in) :: block !< block for computation
      real (kind=RKIND), intent(in) :: minHeight

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

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: statePool, meshPool, tendPool
      real (kind=RKIND), dimension(:), pointer :: sshSubcycleNew
      integer, pointer :: nCellsSolve
      integer :: iCell, k
      integer :: debugUnit
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessTend
      real (kind=RKIND) :: minThickness, layerThick
      character (len=StrKIND) :: debugFilename

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing coefficients
      !
      !-----------------------------------------------------------------

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, timeLevel=1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, timeLevel=2)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
      call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)

      err = 0

      call mpas_log_write( 'Verifying that cells are not dry...')

      ! check to make sure that there is no layer that is too dry
      minThickness = +1.0E34
      do iCell = 1, nCellsSolve
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          ! use ssh as a proxy too for baroclinic mode
          ! Note: wetting-drying currently not supported for either of these time integration methods
          if (trim(config_time_integrator) == 'split_explicit' .or. &
              trim(config_time_integrator) == 'split_implicit' .or. &
              trim(config_time_integrator) == 'split_explicit_ab2') then
            layerThick = min(layerThicknessNew(k, iCell), &
                             (sshSubcycleNew(iCell) + bottomDepth(iCell))/maxLevelCell(iCell))
          else
            layerThick = layerThicknessNew(k, iCell)
          end if
          minThickness = min(minThickness, layerThick)
          if ( layerThick <  minHeight ) then
            ! report error if layerThickness constraint is violated
            if (err == 0) then
              debugFilename = ocn_build_log_filename('mpas_ocean_layerThickness_block_stats_', block % blockID)
              call mpas_new_unit(debugUnit)
              open(unit=debugUnit, file=debugFilename, form='formatted', status='unknown')
            end if
            err = iCell
            write(debugUnit, '(A, I5, A, I5, A, ES14.7, A, ES14.7, A, ES14.7, A)') &
              'ERROR: layerThickness too small at iCell=', iCell, ' k = ', k, &
              ', lon= ', lonCell(iCell), ', lat= ', latCell(iCell), &
              ' with thickness of ', layerThick , '.'
          end if
        end do
      end do

      call mpas_log_write('Minimum thickness is $r.', realArgs=(/minThickness/))

      if ( err > 0) then
        ! end the simulation if layerThickness constraint is violated
        call mpas_log_write( 'Done, some cells have dried.')
        flush(debugUnit)
        close(debugUnit)
        call mpas_release_unit(debugUnit)
        call mpas_log_write( 'ERROR: Layer thickness smaller than $r, see ' // debugFilename , &
                             realArgs=(/ minHeight /), messageType=MPAS_LOG_CRIT)
      else
        call mpas_log_write( 'Done verifying that cells are wet.')
      end if

   end subroutine ocn_wetting_drying_verify !}}}


!***********************************************************************
!
!  routine ocn_prevent_drying_rk4
!
!> \brief   Prevents velocity tendency from causing cell drying
!> \author  Phillip J. Wolfram
!> \date    03/20/2018
!> \details
!>  This routine modifies velocity tendency to prevent cell drying.
!
!-----------------------------------------------------------------------

   subroutine ocn_prevent_drying_rk4(domain, block, dt, rkSubstepWeight, config_zero_drying_velocity, err) !{{{

     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     type (block_type), intent(in) :: block
     real (kind=RKIND), intent(in) :: dt
     real (kind=RKIND), intent(in) :: rkSubstepWeight
     logical, pointer :: config_zero_drying_velocity

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

     integer, intent(out) :: err !< Output: error flag

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------


      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: provisStatePool
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessProvis
      real (kind=RKIND), dimension(:, :), pointer :: normalVelocity

      integer :: iEdge, k

      err = 0

     call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
     call mpas_pool_get_subpool(block % structs, 'state', statePool)
     call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

     call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
     ! use thickness at n because constraint is h_n + dt*T_h > h_min
     call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
     call mpas_pool_get_array(provisStatePool, 'layerThickness', layerThicknessProvis, 1)


     !$omp parallel
     !$omp do schedule(runtime)
     do iEdge = 1, nEdgesAll
       wettingVelocityFactor(:, iEdge) = 0.0_RKIND
     end do
     !$omp end do
     !$omp end parallel

     ! ensure cells stay wet by selectively damping cells with a damping tendency to make
     ! sure tendency doesn't dry cells

     call ocn_wetting_drying_wettingVelocity(domain, layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &
                                             normalTransportVelocity, rkSubstepWeight, wettingVelocityFactor, err)

     ! prevent drying from happening with selective wettingVelocityFactor
     if (config_zero_drying_velocity) then
       !$omp parallel
       !$omp do schedule(runtime) private(k)
       do iEdge = 1, nEdgesAll
         do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)

           if (abs(wettingVelocityFactor(k, iEdge)) > 0.0_RKIND) then
             normalTransportVelocity(k, iEdge) = (1.0_RKIND - &
               wettingVelocityFactor(k, iEdge)) * normalTransportVelocity(k, iEdge)
             normalVelocity(k, iEdge) = (1.0_RKIND - &
               wettingVelocityFactor(k, iEdge)) * normalVelocity(k, iEdge)
           end if

         end do
       end do
       !$omp end do
       !$omp end parallel
     end if

   end subroutine ocn_prevent_drying_rk4 !}}}

!***********************************************************************
!
!  routine ocn_wetting_drying_barotropic_factor
!
!> \brief   Computes velocity to prevent cell drying
!> \author  Phillip J. Wolfram
!> \date    03/19/2018
!> \details
!>  This routine adds wetting velocity opposed to drying motion
!>  to prevent cells from drying.
!
!-----------------------------------------------------------------------
   subroutine ocn_wetting_drying_barotropic_factor(sshNew, btrVelocity, &
                     wettingVelocityBarotropic, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         sshNew, &   !< predicted SSH
         btrVelocity !< barotropic velocity on edges

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(inout) :: &
         wettingVelocityBarotropic !< Input/Output: velocity wettingVelocityFactor

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: cell1, cell2, iEdge, iCell, i

      real (kind=RKIND) :: divFlux, divOutFlux
      real (kind=RKIND) :: wct
      real (kind=RKIND) :: hCrit, hRampMin, hRampMax

      character (len=100) :: log_string

      err = 0

      call mpas_log_write('Barotropic wetdry start')
      wettingVelocityBarotropic(:) = 1.0_RKIND

      if (.not. config_use_wetting_drying .or. &
          .not. config_zero_drying_velocity) then
         call mpas_log_write('Barotropic wetdry return')
         return
      endif

      ! need predicted transport velocity to limit drying flux
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, wct)
      do iCell = 1, nCellsAll
        wct = sshNew(iCell) + bottomDepth(iCell)
        ! if layer thickness is too small, limit divergence flux outwards with
        ! opposite velocity
        hCrit = config_drying_min_cell_height * (maxLevelCell(iCell) - minLevelCell(iCell) + 1) + &
                config_drying_safety_height
        hRampMax = config_zero_drying_velocity_ramp_hmax * &
                   (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
        if (wct <= hCrit) then
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
              wettingVelocityBarotropic(iEdge) = 0.0_RKIND
            end if
          end do
        elseif (config_zero_drying_velocity_ramp .and. &
               (wct > hCrit) .and. &
               (wct <= hRampMax)) then

          hRampMin = config_zero_drying_velocity_ramp_hmin * (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
          ! Following O'Dea et al. (2021), if total upwinded wct is less than
          ! 2*critical thickness, apply damping at each edge
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
              wettingVelocityBarotropic(iEdge) = &
                tanh(50.0_RKIND * (wct - hRampMin)/hRampMin)
            end if
          end do

        end if

      end do
      !$omp end do
      !$omp end parallel
      call mpas_log_write('Barotropic wetdry finish')

   end subroutine ocn_wetting_drying_barotropic_factor !}}}


!***********************************************************************
!
!  routine ocn_wetting_drying_wettingVelocity
!
!> \brief   Computes velocity to prevent cell drying
!> \author  Phillip J. Wolfram
!> \date    03/19/2018
!> \details
!>  This routine adds wetting velocity opposed to drying motion
!>  to prevent cells from drying.
!
!-----------------------------------------------------------------------
   subroutine ocn_wetting_drying_outflux_barotropic_factor(sshCur, btrVelocity, wctEdge, dt, &
                     wettingVelocityBarotropic, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         sshCur, &      !< current SSH
         btrVelocity, & !< barotropic velocity on edges
         wctEdge        !< column thickness on edges

      real (kind=RKIND), intent(in) :: dt !< time step in s

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(inout) :: &
         wettingVelocityBarotropic !< Input/Output: velocity wettingVelocityFactor

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: cell1, cell2, iEdge, iCell, i, cellNei

      real (kind=RKIND) :: divFlux, divOutFlux
      real (kind=RKIND) :: wct, wctNei, sshNei
      real (kind=RKIND) :: hCrit, hRampMin, hRampMax, hRampFactor
      real (kind=RKIND), dimension(:), allocatable :: sshNew

      character (len=100) :: log_string

      err = 0

      call mpas_log_write('Barotropic wetdry start')
      wettingVelocityBarotropic(:) = 1.0_RKIND

      if (.not. config_use_wetting_drying .or. &
          .not. config_zero_drying_velocity) then
         call mpas_log_write('Barotropic wetdry return')
         return
      endif

      allocate(sshNew(nCellsAll))
      do iCell = 1, nCellsAll
        divOutFlux = 0.0_RKIND
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          ! only consider divergence flux leaving the cell
          if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
            divOutFlux = divOutFlux + &
                         btrVelocity(iEdge) * edgeSignOnCell(i, iCell) * &
                         wctEdge(iEdge) * dvEdge(iEdge) * invAreaCell(iCell)
          end if
        end do
        sshNew(iCell) = sshCur(iCell) + dt * divOutFlux
      enddo
      ! need predicted transport velocity to limit drying flux
      hRampFactor = config_zero_drying_velocity_ramp_factor
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, wct)
      do iCell = 1, nCellsAll
        ! Note that this is in contrast to how the baroclinic factor is determined, via only the outgoing flux
        ! Instead this uses the predicted ssh as usual
        wct = sshNew(iCell) + bottomDepth(iCell)
        ! if layer thickness is too small, limit divergence flux outwards with
        ! opposite velocity
        hCrit = (config_drying_min_cell_height + config_drying_safety_height) * &
                (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
                
        hRampMax = config_zero_drying_velocity_ramp_hmax * &
                   (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
        if (wct <= hCrit) then
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
              wettingVelocityBarotropic(iEdge) = 0.0_RKIND
            end if
          end do
          if (config_zero_drying_velocity_gradSSH) then
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              ! Allow flux across edges that are anticipated to be motionless but have an ssh gradient from a wet cell into a dry cell
              if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) == 0.0_RKIND ) then
                cellNei = merge( cellsOnEdge(2, iEdge), cellsOnEdge(1, iEdge), iCell == cellsOnEdge(1, iEdge) )
                if ( cellNei > nCellsAll ) cellNei = iCell ! check with Steven on how to handle this index
                sshNei = sshNew(cellNei)
                wctNei = sshNew(cellNei) + bottomDepth(cellNei)
                if ( wctNei > hCrit .and. sshNew(iCell) < sshNei ) then
                  wettingVelocityBarotropic(iEdge) = 1.0_RKIND
                else
                  wettingVelocityBarotropic(iEdge) = 0.0_RKIND
                end if     
              end if
            end do
          end if     
        elseif (config_zero_drying_velocity_ramp .and. &
               (wct > hCrit) .and. &
               (wct <= hRampMax)) then

          hRampMin = config_zero_drying_velocity_ramp_hmin * (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
          ! Following O'Dea et al. (2021), if total upwinded wct is less than
          ! 2*critical thickness, apply damping at each edge
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
              wettingVelocityBarotropic(iEdge) = &
                tanh(hRampFactor * (wct - hRampMin)/hRampMin)
            end if
          end do
          if (config_zero_drying_velocity_gradSSH) then
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              ! Then go back through all the edges with 0 flux and disallow flux across
              ! edges that do not have an ssh gradient from a wet cell into a dry cell
              if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) == 0.0_RKIND ) then
                cellNei = merge( cellsOnEdge(2, iEdge), cellsOnEdge(1, iEdge), iCell == cellsOnEdge(1, iEdge) )
                if ( cellNei > nCellsAll ) cellNei = iCell ! check with Steven on how to handle this index
                sshNei = sshNew(cellNei)
                wctNei = sshNew(cellNei) + bottomDepth(cellNei)
                if ( wctNei <= hCrit .and. sshNew(iCell) < sshNei ) then
                  wettingVelocityBarotropic(iEdge) = 0.0_RKIND
                end if     
              end if
            end do
          end if

        end if

      end do
      !$omp end do
      !$omp end parallel
      call mpas_log_write('Barotropic wetdry finish')

   end subroutine ocn_wetting_drying_outflux_barotropic_factor !}}}

!***********************************************************************
!
!  routine ocn_wetting_drying_barotropic_factor
!
!> \brief   Computes velocity to prevent cell drying
!> \author  Phillip J. Wolfram
!> \date    03/19/2018
!> \details
!>  This routine adds wetting velocity opposed to drying motion
!>  to prevent cells from drying.
!
!-----------------------------------------------------------------------
   subroutine ocn_wetting_drying_barotropic_factor(sshNew, btrVelocity, &
                     wettingVelocityBarotropic, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         sshNew, &   !< predicted SSH
         btrVelocity !< barotropic velocity on edges

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(inout) :: &
         wettingVelocityBarotropic !< Input/Output: velocity wettingVelocityFactor

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: cell1, cell2, iEdge, iCell, i, cellNei

      real (kind=RKIND) :: divFlux, divOutFlux
      real (kind=RKIND) :: wct, wctNei, sshNei
      real (kind=RKIND) :: hCrit, hRampMin, hRampMax, hRampFactor

      character (len=100) :: log_string

      err = 0

      call mpas_log_write('Barotropic wetdry start')
      wettingVelocityBarotropic(:) = 1.0_RKIND

      if (.not. config_use_wetting_drying .or. &
          .not. config_zero_drying_velocity) then
         call mpas_log_write('Barotropic wetdry return')
         return
      endif

      ! need predicted transport velocity to limit drying flux
      hRampFactor = config_zero_drying_velocity_ramp_factor
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, wct)
      do iCell = 1, nCellsAll
        ! Note that this is in contrast to how the baroclinic factor is determined, via only the outgoing flux
        ! Instead this uses the predicted ssh as usual
        wct = sshNew(iCell) + bottomDepth(iCell)
        ! if layer thickness is too small, limit divergence flux outwards with
        ! opposite velocity
        hCrit = (config_drying_min_cell_height + config_drying_safety_height) * &
                (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
                
        hRampMax = config_zero_drying_velocity_ramp_hmax * &
                   (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
        if (wct <= hCrit) then
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
              wettingVelocityBarotropic(iEdge) = 0.0_RKIND
            end if
          end do
          ! Allow flux across edges that are anticipated to be motionless but have an ssh gradient from a wet cell into a dry cell
          if (config_zero_drying_velocity_gradSSH) then
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) == 0.0_RKIND ) then
                cellNei = merge( cellsOnEdge(2, iEdge), cellsOnEdge(1, iEdge), iCell == cellsOnEdge(1, iEdge) )
                if ( cellNei > nCellsAll ) cellNei = iCell ! check with Steven on how to handle this index
                sshNei = sshNew(cellNei)
                wctNei = sshNew(cellNei) + bottomDepth(cellNei)
                if ( wctNei > hCrit .and. sshNew(iCell) < sshNei ) then
                  wettingVelocityBarotropic(iEdge) = 1.0_RKIND
                else
                  wettingVelocityBarotropic(iEdge) = 0.0_RKIND
                end if     
              end if
            end do
          end if
        elseif (config_zero_drying_velocity_ramp .and. &
               (wct > hCrit) .and. &
               (wct <= hRampMax)) then

          hRampMin = config_zero_drying_velocity_ramp_hmin * (maxLevelCell(iCell) - minLevelCell(iCell) + 1)
          ! Following O'Dea et al. (2021), if total upwinded wct is less than
          ! 2*critical thickness, apply damping at each edge
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
              wettingVelocityBarotropic(iEdge) = &
                tanh(hRampFactor * (wct - hRampMin)/hRampMin)
            end if
          end do
          if (config_zero_drying_velocity_gradSSH) then
            ! Then go back through all the edges with 0 flux and disallow flux across
            ! edges that do not have an ssh gradient from a wet cell into a dry cell
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if ( btrVelocity(iEdge) * edgeSignOnCell(i, iCell) == 0.0_RKIND ) then
                cellNei = merge( cellsOnEdge(2, iEdge), cellsOnEdge(1, iEdge), iCell == cellsOnEdge(1, iEdge) )
                if ( cellNei > nCellsAll ) cellNei = iCell ! check with Steven on how to handle this index
                sshNei = sshNew(cellNei)
                wctNei = sshNew(cellNei) + bottomDepth(cellNei)
                if ( wctNei <= hCrit .and. sshNew(iCell) < sshNei ) then
                  wettingVelocityBarotropic(iEdge) = 0.0_RKIND
                end if     
              end if
            end do
          end if

        end if

      end do
      !$omp end do
      !$omp end parallel
      call mpas_log_write('Barotropic wetdry finish')

   end subroutine ocn_wetting_drying_barotropic_factor !}}}

!***********************************************************************
!
!  routine ocn_wetting_drying_full_depth_factor
!
!> \brief   Computes velocity to prevent cell drying
!> \author  Phillip J. Wolfram
!> \date    03/19/2018
!> \details
!>  This routine adds wetting velocity opposed to drying motion
!>  to prevent cells from drying.
!
!-----------------------------------------------------------------------
   subroutine ocn_wetting_drying_full_depth_factor(layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &

       normalVelocity, dt, wettingVelocityFactor, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessCur    !< Input: layer thickness at old time

       real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessProvis    !< Input: provisional layer thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeFlux  !< Input: flux-related layerThickness at an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport

      real (kind=RKIND), intent(in) :: &
         dt     !< Input: time step

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         wettingVelocityFactor          !< Input/Output: velocity wettingVelocityFactor

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: cell1, cell2, cellNei, iEdge, iCell, k, i

      real (kind=RKIND) :: divFlux, divOutFlux
      real (kind=RKIND) :: layerThickness
      real (kind=RKIND) :: hCrit, hRampMin

      character (len=100) :: log_string

      err = 0

      if (.not. config_zero_drying_velocity) return

      allocate(columnThickness(nCellsAll))
      ! need predicted transport velocity to limit drying flux
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, k, divOutFlux, layerThickness)
      do iCell = 1, nCellsAll
        columnThickness(iCell) = 0.0_RKIND
        hCrit = 0.0_RKIND
        hRampMax = 0.0_RKIND
        hRampMin = 0.0_RKIND

        do k = minLevelCell(iCell), maxLevelCell(iCell)
          divOutFlux = 0.0_RKIND

          layerThickness = min(layerThicknessProvis(k, iCell), layerThicknessCur(k, iCell))
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeTop(iEdge) .and. k >= minLevelEdgeBot(iEdge)) then
              ! only consider divergence flux leaving the cell
              if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
                divOutFlux = divOutFlux + &
                             normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) * &
                             layerThickEdgeFlux(k, iEdge) * dvEdge(iEdge) * &
                             invAreaCell(iCell)
              end if
            end if
          end do
          layerThickness = layerThickness + dt * divOutFlux
          columnThickness(iCell) = columnThickness(iCell) + layerThickness
          hCrit = hCrit + config_drying_min_cell_height + config_drying_safety_height
          hRampMax = hRampMax + config_zero_drying_velocity_ramp_hmax
          hRampMin = hRampMin + config_zero_drying_velocity_ramp_hmin
        enddo
        ! The projected column thickness only accounts for outgoing fluxes
        ! if column thickness is too small, limit flux throughout the full column
        if (columnThickness(iCell) <= hCrit) then
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND
                end if
              end if
            end do
          end do
        end if
        if (config_zero_drying_velocity_gradSSH) then
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
                ! Allow flux across edges that are anticipated to be motionless but have an ssh gradient from a wet cell into a dry cell
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) == 0.0_RKIND ) then
                  cellNei = merge( cellsOnEdge(2, iEdge), cellsOnEdge(1, iEdge), iCell == cellsOnEdge(1, iEdge) )
                  if ( cellNei > nCellsAll ) cellNei = iCell ! check with Steven on how to handle this index
                  sshCur = columnThickness(iCell) - bottomDepth(iCell)
                  sshNei = columnThickness(cellNei) - bottomDepth(cellNei)
                  if ( columnThickness(cellNei) > hCrit .and. sshCur < sshNei ) then
                    wettingVelocityFactor(k, iEdge) = 0.0_RKIND
                  else
                    wettingVelocityFactor(k, iEdge) = 1.0_RKIND
                  end if     
                end if
              end if
            end do
          end do
        end if
        if (config_zero_drying_velocity_ramp) then
          hRampFactor = config_zero_drying_velocity_ramp_factor
          if ((columnThickness(iCell) > hCrit ) .and. &
              (columnThickness(iCell) <= hRampMax)) then

            ! Following O'Dea et al. (2021), if total upwinded wct is less than
            ! 2*critical thickness, apply damping at each edge
            do k = minLevelCell(iCell), maxLevelCell(iCell)
              do i = 1, nEdgesOnCell(iCell)
                iEdge = edgesOnCell(i, iCell)
                if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
                  if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                    wettingVelocityFactor(k, iEdge) = 1.0_RKIND - &
                      tanh(hRampFactor * (columnThickness(iCell) - hRampMin)/hRampMin)
                  end if
                end if
              end do
            end do
            if (config_zero_drying_velocity_gradSSH) then
              do k = minLevelCell(iCell), maxLevelCell(iCell)
                do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  ! Then go back through all the edges with 0 flux and disallow flux when the neighboring cell is dry
                  ! but the ssh gradient would favor flow from the neighboring dry cell
                  if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) == 0.0_RKIND ) then
                    cellNei = merge( cellsOnEdge(2, iEdge), cellsOnEdge(1, iEdge), iCell == cellsOnEdge(1, iEdge) )
                    if ( cellNei > nCellsAll ) cellNei = iCell ! check with Steven on how to handle this index
                    sshCur = columnThickness(iCell) - bottomDepth(iCell)
                    sshNei = columnThickness(cellNei) - bottomDepth(cellNei)
                    if ( columnThickness(cellNei) <= hCrit .and. sshCur < sshNei ) then
                      wettingVelocityFactor(k, iEdge) = 1.0_RKIND
                    end if     
                  end if
                end do
              end do
            end if
          end if
        end if

      end do
      !$omp end do
      !$omp end parallel

   end subroutine ocn_wetting_drying_full_depth_factor !}}}

!***********************************************************************
!
!  routine ocn_wetting_drying_wettingVelocity
!
!> \brief   Computes velocity to prevent cell drying
!> \author  Phillip J. Wolfram
!> \date    03/19/2018
!> \details
!>  This routine adds wetting velocity opposed to drying motion
!>  to prevent cells from drying.
!
!-----------------------------------------------------------------------

   subroutine ocn_wetting_drying_wettingVelocity(domain, layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &
                                                 normalVelocity, dt, wettingVelocityFactor, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessCur    !< Input: layer thickness at old time

       real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessProvis    !< Input: provisional layer thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeFlux  !< Input: flux-related layerThickness at an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport

      real (kind=RKIND), intent(in) :: &
         dt     !< Input: time step

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), optional, intent(inout) :: &
         domain !< Input/Output: domain information, needed for halo exchange


      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         wettingVelocityFactor          !< Input/Output: velocity wettingVelocityFactor

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: cell1, cell2, iEdge, iCell, k, i, cellNei

      real (kind=RKIND) :: divFlux, divOutFlux
      real (kind=RKIND) :: sshCur, sshNei
      real (kind=RKIND) :: hCrit, hRampMin, hRampFactor
      real (kind=RKIND), dimension(:), allocatable :: columnThickness
      real (kind=RKIND), dimension(:,:), allocatable :: layerThickness
      character (len=100) :: log_string

      integer:: cellDummy(2), cellCur, CellNei
      real (kind=RKIND) :: sshCur, sshNei
      real (kind=RKIND), dimension(:), pointer :: sshCell
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCell


      err = 0

      if (.not. config_zero_drying_velocity) return

      allocate(columnThickness(nCellsAll))
      allocate(layerThickness(nVertLevels, nCellsAll))
      ! need predicted transport velocity to limit drying flux
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, k, divOutFlux)
      do iCell = 1, nCellsAll
        columnThickness(iCell) = 0.0_RKIND
        layerThickness(:, iCell) = 0.0_RKIND
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          divOutFlux = 0.0_RKIND
          layerThickness(k, iCell) = min(layerThicknessProvis(k, iCell), layerThicknessCur(k, iCell))
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeTop(iEdge) .and. k >= minLevelEdgeBot(iEdge)) then
              ! only consider divergence flux leaving the cell
              if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
                divOutFlux = divOutFlux + &
                             normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) * &
                             layerThickEdgeFlux(k, iEdge) * dvEdge(iEdge) * &
                             invAreaCell(iCell)
              end if
            end if
          end do
          layerThickness(k, iCell) = layerThickness(k, iCell) + dt * divOutFlux
          columnThickness(iCell) = columnThickness(iCell) + layerThickness(k, iCell)

          ! if layer thickness is too small, limit divergence flux outwards with
          ! opposite velocity
          if (layerThickness(k, iCell) <= &
              hCrit + config_drying_safety_height) then
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND
                end if
              end if
            end do
          elseif (config_zero_drying_velocity_ramp .and. &
                 (layerThickness(k, iCell) > &
                 hCrit + config_drying_safety_height) .and. &
                 (layerThickness(k, iCell) <= config_zero_drying_velocity_ramp_hmax)) then

            hRampMin = config_zero_drying_velocity_ramp_hmin
            hRampFactor = config_zero_drying_velocity_ramp_factor
            ! Following O'Dea et al. (2021), if total upwinded wct is less than
            ! 2*critical thickness, apply damping at each edge
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND - &
                    tanh(hRampFactor * (layerThickness(k, iCell) - hRampMin)/hRampMin)
                end if
              end if
            end do

          layerThicknessCell(k, iCell) = min(layerThicknessProvis(k, iCell), layerThicknessCur(k, iCell)) + dt * divOutFlux
          if ( .not. config_use_ssh_gradient_wetting_drying ) then
             call ocn_wetting_velocity_factor_on_cell_edges( &
                      wettingVelocityFactor, layerThicknessCell(k, iCell), normalVelocity, iCell, k)
          end if
        end do
      end do
      !$omp end do
      !$omp end parallel

      if (config_zero_drying_velocity_gradSSH) then
        !$omp parallel
        !$omp do schedule(runtime) private(i, iEdge, k)
        do iCell = 1, nCellsAll
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeTop(iEdge) .and. k >= minLevelEdgeBot(iEdge) .and. &
                  wettingVelocityFactor(k, iEdge) .ne. 0.0_RKIND .and. &
                  normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) == 0.0_RKIND ) then
                ! Then go back through all the edges with 0 flux and disallow flux when the neighboring cell is dry
                ! but the ssh gradient would favor flow from the neighboring dry cell
                cellNei = merge( cellsOnEdge(2, iEdge), cellsOnEdge(1, iEdge), iCell == cellsOnEdge(1, iEdge) )
                if ( cellNei > nCellsAll ) cellNei = iCell ! check with Steven on how to handle this index
                sshCur = columnThickness(iCell) - bottomDepth(iCell)
                sshNei = columnThickness(cellNei) - bottomDepth(cellNei)
                if (config_zero_drying_velocity_ramp) then
                  if ( layerThickness(k, cellNei) <= hCrit .and. sshCur < sshNei ) then
                    wettingVelocityFactor(k, iEdge) = 1.0_RKIND
                  end if
                else
                  if ( layerThickness(k, cellNei) > hCrit .and. sshCur < sshNei ) then
                    wettingVelocityFactor(k, iEdge) = 0.0_RKIND
                  else
                    wettingVelocityFactor(k, iEdge) = 1.0_RKIND
                  end if
                end if
              end if
            end do
          end do
        end do
        !$omp end do
        !$omp end parallel
      end if

   end subroutine ocn_wetting_drying_wettingVelocity !}}}

!***********************************************************************
!
!  routine ocn_wetting_velocity_factor_on_cell_edges
!
!> \brief   Computes wettingVelocityFactor at edges of iCell
!> \author  Giacomo Capodaglio
!> \date    09/06/2023
!> \details
!>  This routine computes wettingVelocityFactor at the edges of iCell
!
!-----------------------------------------------------------------------
   subroutine ocn_wetting_velocity_factor_on_cell_edges(wettingVelocityFactor, layerThicknessTmp, normalVelocity, iCell, k)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         layerThicknessTmp    !< Input: layer thickness to use for limiting

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity       !< Input: transport

      integer, intent(in) :: &
         iCell                !< Input: wettingVelocityFactor is computed at the edges of iCell

      integer, intent(in) :: &
         k                    !< Input: wettingVelocityFactor is computed at layer k

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         wettingVelocityFactor          !< Input/Output: velocity wettingVelocityFactor

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i, iEdge

      real (kind=RKIND) :: hCrit, hRampMin

      hCrit = config_drying_min_cell_height

      ! if layer thickness is too small, limit divergence flux outwards with
      ! opposite velocity
      if (layerThicknessTmp <= &
         hCrit + config_drying_safety_height) then
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
               if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND
               end if
            end if
         end do
      elseif (config_zero_drying_velocity_ramp .and. &
              (layerThicknessTmp > &
              hCrit + config_drying_safety_height) .and. &
              (layerThicknessTmp <= config_zero_drying_velocity_ramp_hmax)) then

         hRampMin = config_zero_drying_velocity_ramp_hmin
         ! Following O'Dea et al. (2021), if total upwinded wct is less than
         ! 2*critical thickness, apply damping at each edge
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
               if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND - &
                  tanh(50.0_RKIND * (layerThicknessTmp - hRampMin)/hRampMin)
               end if
            end if
         end do

      end if

   end subroutine ocn_wetting_velocity_factor_on_cell_edges!}}}

end module ocn_wetting_drying

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
