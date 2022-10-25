!==============================================================================
! MODULE: utilities_module
!
! Authors
! -------
! Dr. Randal J. Barnes
! Department of Civil, Environmental, and Geo- Engineering
! University of Minnesota
! Minneapolis, MN, 55455
!
! Dr. Elizabeth Barnes
! Department of Atmospheric Science
! Colorado State University
! Fort Collins, CO, 80525
!
! Dr. Mark DeMaria
! Cooperative Institute for Research in the Atmosphere
! Colorado State University
! Fort Collins, CO, 80523
!
! Version
! -------
! * 25 October 2022
!
!==============================================================================
module utilities_module
   implicit none

   !----------------------------------------------------------------------------
   ! Visibility
   !----------------------------------------------------------------------------
   private

   public :: to_upper
   public :: to_lower
   public :: get_time_string
   public :: get_date_string

contains

   !-----------------------------------------------------------------------------
   ! FUNCTION: to_upper
   !
   !  Convert all alphabetic characters in a string to upper case.
   !
   ! Arguments
   ! ---------
   ! string_in : character(len=*)
   !
   ! Returns
   ! -------
   ! string_out : character(len=len(string_in))
   !
   !-----------------------------------------------------------------------------
   function to_upper(string_in) result(string_out)
      character(len=*), intent(in)  :: string_in
      character(len=len(string_in)) :: string_out
      integer :: i, j

      ! Shift all alphabetic characters to upper case.
      do i = 1, len(string_in)
         j = iachar(string_in(i:i))
         if (j >= iachar("a") .and. j <= iachar("z")) then
            string_out(i:i) = achar(j - 32)
         else
            string_out(i:i) = string_in(i:i)
         end if
      end do

   end function to_upper

   !-----------------------------------------------------------------------------
   ! FUNCTION: to_lower
   !
   !  Convert all alphabetic characters in a string to lower case.
   !
   ! Arguments
   ! ---------
   ! string_in : character(len=*)
   !
   ! Returns
   ! -------
   ! string_out : character(len=len(string_in))
   !
   !-----------------------------------------------------------------------------
   function to_lower(string_in) result(string_out)
      character(len=*), intent(in)  :: string_in
      character(len=len(string_in)) :: string_out
      integer :: i, j

      ! Shift all alphabetic characters to upper case.
      do i = 1, len(string_in)
         j = iachar(string_in(i:i))
         if (j >= iachar("A") .and. j <= iachar("Z")) then
            string_out(i:i) = achar(j + 32)
         else
            string_out(i:i) = string_in(i:i)
         end if
      end do

   end function to_lower

   !-----------------------------------------------------------------------------
   ! FUNCTION: get_time_string
   !
   !  Returns a formatted time string.
   !-----------------------------------------------------------------------------
   function get_time_string() result(buffer)
      character(12) :: buffer

      character(10) :: system_time

      call date_and_time(time = system_time)
      write(buffer, "('[', A2, ':', A2, ':', A4, ']')") &
         system_time(1:2), system_time(3:4), system_time(5:8)

   end function get_time_string

   !-----------------------------------------------------------------------------
   ! FUNCTION: get_date_string
   !
   !  Returns a formatted date string.
   !-----------------------------------------------------------------------------
   function get_date_string() result(buffer)
      character(12) :: Buffer

      character(8) :: system_date

      call date_and_time(date = system_date)
      write(buffer, "('[', A4, ':', A2, ':', A2, ']')") &
         system_date(1:4), system_date(5:6), system_date(7:8)

   end function get_date_string

end module utilities_module
