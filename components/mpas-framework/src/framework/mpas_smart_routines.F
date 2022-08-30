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

   public :: mpas_smart_init

#  include "enum_fortran.inc"

   contains

   subroutine mpas_smart_init(client, ierr)

      implicit none

      type(client_type), intent(out) :: client
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
            write(log_string,*) 'FWK: Error initializing redis client:',err_client
            call mpas_log_write(log_string)
            ierr = ior(ierr, 1)
         else
            call mpas_log_write('FWK: Client initialized')
         endif
      endif

   end subroutine mpas_smart_init

   !subroutine mpas_smart_stream_send(client, stream, ierr)

   !   implicit none

   !   type(client_type), intent(in) :: client
   !   !type (MPAS_streamManager_type), intent(inout) :: manager
   !   !character (len=*), intent(in), optional :: streamID
   !   type (MPAS_Stream_type), intent(in) :: stream
   !   !type (MPAS_stream_list_type), pointer :: stream_cursor
   !   integer, intent(out), optional :: ierr

   !   type (field_list_type), pointer :: field_cursor
   !   type (field1dReal), pointer :: field_1d_ptr
   !   !call MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=ierr))
   !   !
   !   ! Loop over fields in the stream
   !   !
   !   field_cursor => stream % fieldList
   !   do while (associated(field_cursor))
   !      if (field_cursor % field_type == FIELD_1D_REAL) then
   !         field_1d_ptr = field_cursor % real1dField
   !         call mpas_smart_field_send_1d(client, field_1d_ptr, ierr)
   !      else
   !         call mpas_log_write('field is not of type 1D REAL')
   !      endif
   !   end do

   !end subroutine mpas_smart_stream_send

   subroutine mpas_smart_field_send(client, domain, fieldName, ierr)

      type(client_type), intent(in) :: client
      type (domain_type), intent(in) :: domain
      character(len=100), intent(in) :: fieldName
      integer, intent(out), optional :: ierr

      type (field1dReal), pointer :: field_ptr
      character(len=100) :: log_string
      integer(kind=enum_kind) :: err_client

      call mpas_pool_get_field_info(domain%blocklist%allFields, &
                                    smartFieldName, fieldInfo)
      WRITE(log_string,*) smartFieldName, 'has type', fieldInfo%fieldType
      call mpas_log_write(log_string)
      if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
         call mpas_pool_get_field(domain%blocklist%allFields, &
                                  smartFieldName, sendArrayReal1d)
         WRITE(log_string,*) smartFieldName, 'retrieved from stream', &
                             sendArrayReal1d
         call mpas_smart_field_send_1d(client, field_ptr, ierr)
      endif

   end subroutine mpas_smart_field_send

   subroutine mpas_smart_field_send_1d(client, field_ptr, ierr)

      implicit none

      type(client_type), intent(in) :: client
      type (field1dReal), pointer, intent(in) :: field_ptr
      integer, intent(out), optional :: ierr

      character(len=100) :: log_string
      integer(kind=enum_kind) :: err_client
      real (kind=RKIND), dimension(:), allocatable :: field_temp
      integer, dimension(1) :: field_dims

      call mpas_pool_get_dimension(field_ptr % block % dimensions, &
                                   trim(field_ptr % dimNames(1), dim1)
      allocate(field_temp(field_dims(1)))
      field_temp = field_ptr % array

      err_client = client%put_tensor("send_array", field_temp, &
                                     shape(sendArrayReal1d))
      if (err_client /= SRNoError) then
         WRITE(log_string,*) 'Error put tensor:', smartFieldName
         call mpas_log_write(log_string)
         ierr = ior(ierr, 1)
      endif

   end subroutine mpas_smart_field_send
end module mpas_smart_routines
