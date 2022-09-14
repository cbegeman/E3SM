! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_smart_routines

#define COMMA ,
#define STREAM_DEBUG_WRITE(M) call mpas_log_write(M)

#ifdef SINGLE_PRECISION
#define REAL_IO_TYPE MPAS_IO_REAL
#else
#define REAL_IO_TYPE MPAS_IO_DOUBLE
#endif

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_field_routines
   use mpas_threading
   use mpas_log
   use smartredis_client, only: client_type

   implicit none
   save

   public :: mpas_smart_init

#  include "enum_fortran.inc"

! --- Private module variables ---
   type(client_type) :: client

   interface mpas_smart_field_client_send
      module procedure mpas_smart_field_client_send_1d
      module procedure mpas_smart_field_client_send_2d
      module procedure mpas_smart_field_client_send_3d
      module procedure mpas_smart_field_client_send_4d
   end interface

   interface mpas_smart_field_client_recv
      module procedure mpas_smart_field_client_recv_1d
      module procedure mpas_smart_field_client_recv_2d
      module procedure mpas_smart_field_client_recv_3d
      module procedure mpas_smart_field_client_recv_4d
   end interface

   contains

   subroutine mpas_smart_init(client_out, ierr)

      !implicit none

      type(client_type), intent(out) :: client_out
      integer, intent(out), optional :: ierr

      character(len=100) :: log_string
      integer(kind=enum_kind) :: err_client

      call get_environment_variable('SSDB',log_string)
      call mpas_log_write('SSDB='//TRIM(log_string))
      if (len(TRIM(log_string)) < 5) then
         call mpas_log_write('FWK: No SSDB')!, &
      else
         err_client = client%initialize(.false.)
         if (err_client /= SRNoError) then
            write(log_string,*) 'FWK: Error initializing redis client:', &
                                 err_client
            call mpas_log_write(log_string)
            ierr = ior(ierr, 1)
         else
            call mpas_log_write('FWK: Client initialized')
            client_out = client
         endif
      endif

   end subroutine mpas_smart_init

   subroutine mpas_smart_database_init(ierr)

      !implicit none

      integer, intent(out), optional :: ierr

      character(len=100) :: log_string
      integer(kind=enum_kind) :: err_client

      if (err_client /= SRNoError) then
         write(log_string,*) 'FWK: Error initializing redis database:', &
                              err_client
         call mpas_log_write(log_string)
         ierr = ior(ierr, 1)
      else
         call mpas_log_write('FWK: Database initialized')
      endif

   end subroutine mpas_smart_database_init

!   --------------------------------------------------------------------
!    SUBROUTINE mpas_smart_stream_send
!      Send the stream to a client.  
!   --------------------------------------------------------------------
!   subroutine mpas_smart_stream_send(stream, ierr)
 
!      implicit none
 
!      !type (MPAS_streamManager_type), intent(inout) :: manager
!      !character (len=*), intent(in), optional :: streamID
!      type (MPAS_Stream_type), intent(in) :: stream
!      !type (MPAS_stream_list_type), pointer :: stream_cursor
!      integer, intent(out), optional :: ierr
 
!      type (field_list_type), pointer :: field_cursor
!      type (field1dReal), pointer :: field_1d_ptr
!      !call MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=ierr))
!      !
!      ! Loop over fields in the stream
!      !
!      field_cursor => stream % fieldList
!      do while (associated(field_cursor))
!         if (field_cursor % field_type == FIELD_1D_REAL) then
!            field_1d_ptr = field_cursor % real1dField
!            call mpas_smart_field_send_1d(client, field_1d_ptr, ierr)
!         else
!            call mpas_log_write('field is not of type 1D REAL')
!         endif
!      end do
 
!   end subroutine mpas_smart_stream_send
 
!   --------------------------------------------------------------------
!    SUBROUTINE mpas_smart_stream_database_send
!      Send the stream to a database.  
!      TODO: if decomposed and blockWrite, send each block to a 
!      different database
!   --------------------------------------------------------------------

!   subroutine mpas_smart_stream_database_send(stream, database_list)
!   
!      type (field_list_type), pointer :: field_cursor
 
!      field_cursor => stream % fieldList
!      do while (associated(field_cursor))
 
!               if (field_cursor % isDecomposed) then
!                  ! Gather field from across multiple blocks
!                  field_2dreal_ptr => field_cursor % real2dField
!                  i = 1
!                  do while (associated(field_2dreal_ptr))
!                     i = i + ownedSize
!                     if ( .not. stream % blockWrite ) then
!                        field_2dreal_ptr => field_2dreal_ptr % next
!                     else
!                        nullify(field_2dreal_ptr)
!                     end if
!                  end do
!               end if
!               if ( .not. stream % blockWrite ) then
 
!   end subroutine mpas_smart_database_send

   subroutine mpas_smart_config_send(domain, configName, ierr)

      type (domain_type), intent(in) :: domain
      character(len=*), intent(in) :: configName
      integer, intent(out), optional :: ierr

      type (mpas_pool_field_info_type) :: fieldInfo
      real (kind=RKIND), pointer :: config_ptr
      real (kind=RKIND), dimension(1) :: config_temp
      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client

      call mpas_pool_get_config(domain % configs, &
                                trim(configName), &
                                config_ptr)
      config_temp(1) = config_ptr
      err_client = client%put_tensor(trim(configName), &
                                     config_temp, shape(config_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'FWK: Error put config:', trim(configName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'FWK: Put config:', trim(configName)
         call mpas_log_write(trim(log_string))
      endif

   end subroutine mpas_smart_config_send

   subroutine mpas_smart_config_recv(domain, configName, ierr)

      type (domain_type), intent(in) :: domain
      character(len=*), intent(in) :: configName
      integer, intent(out), optional :: ierr

      type (mpas_pool_field_info_type) :: fieldInfo
      real (kind=RKIND), pointer :: config_ptr
      real (kind=RKIND), dimension(1) :: config_temp
      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client

      call mpas_pool_get_config(domain % configs, &
                                trim(configName), &
                                config_ptr)
      err_client = client%unpack_tensor(trim(configName), &
                                        config_temp, shape(config_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'FWK: Error unpack config:', trim(configName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'FWK: unpack config:', trim(configName)
         call mpas_log_write(trim(log_string))
      endif
      config_ptr = config_temp(1)

   end subroutine mpas_smart_config_recv

   subroutine mpas_smart_field_send(domain, fieldName, ierr)

      type (domain_type), intent(in) :: domain
      character(len=StrKIND), intent(in) :: fieldName
      integer, intent(out), optional :: ierr

      type (mpas_pool_field_info_type) :: fieldInfo
      type (field1DReal), pointer :: field_real1d_ptr
      type (field2DReal), pointer :: field_real2d_ptr
      type (field3DReal), pointer :: field_real3d_ptr
      type (field4DReal), pointer :: field_real4d_ptr
      character(len=StrKIND) :: log_string

      call mpas_pool_get_field_info(domain%blocklist%allFields, &
                                    fieldName, fieldInfo)
      call mpas_log_write(log_string)
      if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
         if ( fieldInfo % nDims == 1 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real1d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_send(field_real1d_ptr, ierr)
         elseif ( fieldInfo % nDims == 2 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real2d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_send(field_real2d_ptr, ierr)
         elseif ( fieldInfo % nDims == 3 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real3d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_send(field_real3d_ptr, ierr)
         elseif ( fieldInfo % nDims == 4 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real4d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_send(field_real4d_ptr, ierr)
         else
            call mpas_log_write(trim(fieldName)// &
                                ' is not of a supported dimension')
         endif
      else
         call mpas_log_write(trim(fieldName)//' is not real')
      endif

   end subroutine mpas_smart_field_send

   subroutine mpas_smart_field_recv(domain, fieldName, ierr)

      type (domain_type), intent(in) :: domain
      character(len=StrKIND), intent(in) :: fieldName
      integer, intent(out), optional :: ierr

      type (mpas_pool_field_info_type) :: fieldInfo
      type (field1DReal), pointer :: field_real1d_ptr
      type (field2DReal), pointer :: field_real2d_ptr
      type (field3DReal), pointer :: field_real3d_ptr
      type (field4DReal), pointer :: field_real4d_ptr
      character(len=StrKIND) :: log_string

      call mpas_pool_get_field_info(domain%blocklist%allFields, &
                                    fieldName, fieldInfo)
      call mpas_log_write(log_string)
      if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
         if ( fieldInfo % nDims == 1 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real1d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_recv(field_real1d_ptr, ierr)
         elseif ( fieldInfo % nDims == 2 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real2d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_recv(field_real2d_ptr, ierr)
         elseif ( fieldInfo % nDims == 3 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real3d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_recv(field_real3d_ptr, ierr)
         elseif ( fieldInfo % nDims == 4 ) then
            call mpas_pool_get_field(domain%blocklist%allFields, &
                                     trim(fieldName), field_real4d_ptr)
            call mpas_log_write(trim(fieldName)//' retrieved')
            call mpas_smart_field_client_recv(field_real4d_ptr, ierr)
         else
            call mpas_log_write(trim(fieldName)// &
                                ' is not of a supported dimension')
         endif
      else
         call mpas_log_write(trim(fieldName)//' is not real')
      endif

   end subroutine mpas_smart_field_recv


!  INTERFACE mpas_smart_field_client_send
!  TODO add optional argument to provide a different key
   subroutine mpas_smart_field_client_send_1d(field_ptr, ierr)

      type (field1dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:), allocatable :: field_temp
      integer, dimension(1) :: field_dims

      field_dims = field_ptr % dimSizes
      !call mpas_pool_get_dimension(field_ptr % block % dimensions, &
      !                             trim(field_ptr % dimNames(1)), field_dims(1))
      allocate(field_temp(field_dims(1)))
      field_temp = field_ptr % array

      err_client = client%put_tensor(trim(field_ptr % fieldName), &
                                     field_temp, &
                                     shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'FWK: Error put tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'FWK: Put tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_send_1d

   subroutine mpas_smart_field_client_send_2d(field_ptr, ierr)

      type (field2dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:,:), allocatable :: field_temp
      integer, dimension(2) :: field_dims

      field_dims = field_ptr % dimSizes
      allocate(field_temp(field_dims(1), field_dims(2)))
      field_temp = field_ptr % array

      err_client = client%put_tensor(trim(field_ptr % fieldName), &
                                     field_temp, &
                                     shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'Error put tensor:', &
                             trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'Put tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_send_2d

   subroutine mpas_smart_field_client_send_3d(field_ptr, ierr)

      type (field3dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:,:,:), allocatable :: field_temp
      integer, dimension(3) :: field_dims

      field_dims = field_ptr % dimSizes
      allocate(field_temp(field_dims(1), field_dims(2), field_dims(3)))
      field_temp = field_ptr % array

      err_client = client%put_tensor(trim(field_ptr % fieldName), &
                                     field_temp, &
                                     shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'Error put tensor:', &
                             trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'Put tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_send_3d

   subroutine mpas_smart_field_client_send_4d(field_ptr, ierr)

      type (field4dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:,:,:,:), allocatable :: field_temp
      integer, dimension(4) :: field_dims

      field_dims = field_ptr % dimSizes
      allocate(field_temp(field_dims(1), field_dims(2), field_dims(3), &
               field_dims(4)))
      field_temp = field_ptr % array

      err_client = client%put_tensor(trim(field_ptr % fieldName), &
                                     field_temp, &
                                     shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'Error put tensor:', &
                             trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'Put tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_send_4d

!  INTERFACE mpas_smart_field_client_recv
!  TODO add optional argument to provide a different key
   subroutine mpas_smart_field_client_recv_1d(field_ptr, ierr)

      type (field1dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:), allocatable :: field_temp
      integer, dimension(1) :: field_dims

      field_dims = field_ptr % dimSizes
      allocate(field_temp(field_dims(1)))

      err_client = client%unpack_tensor(trim(field_ptr % fieldName), &
                                     field_temp, &
                                     shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'FWK: Error unpack tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'FWK: Unpack tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      field_ptr % array = field_temp
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_recv_1d

   subroutine mpas_smart_field_client_recv_2d(field_ptr, ierr)

      type (field2dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:,:), allocatable :: field_temp
      integer, dimension(2) :: field_dims

      field_dims = field_ptr % dimSizes
      allocate(field_temp(field_dims(1), field_dims(2)))

      err_client = client%unpack_tensor(trim(field_ptr % fieldName), &
                                        field_temp, &
                                        shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'Error unpack tensor:', &
                             trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'FWK: unpack tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      field_ptr % array = field_temp
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_recv_2d

   subroutine mpas_smart_field_client_recv_3d(field_ptr, ierr)

      type (field3dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:,:,:), allocatable :: field_temp
      integer, dimension(3) :: field_dims

      field_dims = field_ptr % dimSizes
      allocate(field_temp(field_dims(1), field_dims(2), field_dims(3)))

      err_client = client%unpack_tensor(trim(field_ptr % fieldName), &
                                        field_temp, &
                                        shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'Error unpack tensor:', &
                             trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'FWK: unpack tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      field_ptr % array = field_temp
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_recv_3d

   subroutine mpas_smart_field_client_recv_4d(field_ptr, ierr)

      type (field4dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=StrKIND) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:,:,:,:), allocatable :: field_temp
      integer, dimension(4) :: field_dims

      field_dims = field_ptr % dimSizes
      allocate(field_temp(field_dims(1), field_dims(2), field_dims(3), &
               field_dims(4)))

      err_client = client%unpack_tensor(trim(field_ptr % fieldName), &
                                        field_temp, &
                                        shape(field_temp))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'Error unpack tensor:', &
                             trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
         ierr = ior(ierr, 1)
      else
         WRITE(log_string,*) 'Unpack tensor:', trim(field_ptr % fieldName)
         call mpas_log_write(trim(log_string))
      endif
      field_ptr % array = field_temp
      deallocate(field_temp)

   end subroutine mpas_smart_field_client_recv_4d

end module mpas_smart_routines
