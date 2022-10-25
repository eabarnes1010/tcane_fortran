!==============================================================================
! MODULE: isclose_module
!
! Code that supplies the Fortran 2018 "isclose" functionality to pre-2018
! compilers.
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
module isclose_module
   implicit none

   !----------------------------------------------------------------------------
   ! Visibility
   !----------------------------------------------------------------------------
   private

   public :: isclose
   public :: isclose_matrix
   public :: isclose_scalar
   public :: isclose_vector

   !-----------------------------------
   ! Interfaces
   !-----------------------------------
   interface isclose
      procedure :: isclose_matrix
      procedure :: isclose_scalar
      procedure :: isclose_vector
   end interface

contains

   !---------------------------------------------------------------------------
   ! isclose_matrix
   !---------------------------------------------------------------------------
   function isclose_matrix(actual, desired, rtol, atol) result( flag )
      real(8), dimension(:,:), intent(in)           :: actual
      real(8), dimension(:,:), intent(in)           :: desired
      real(8),                 intent(in), optional :: rtol
      real(8),                 intent(in), optional :: atol
      logical                                       :: flag

      flag = .true.

      if (present(rtol)) then
         flag = flag .and. maxval(abs(actual - desired)/abs(desired)) <= rtol
      end if

      if (present(atol)) then
         flag = flag .and. maxval(abs(actual - desired)) <= atol
      end if
   end function isclose_matrix

   !---------------------------------------------------------------------------
   ! isclose_scalar
   !---------------------------------------------------------------------------
   function isclose_scalar(actual, desired, rtol, atol) result( flag )
      real(8), intent(in)           :: actual
      real(8), intent(in)           :: desired
      real(8), intent(in), optional :: rtol
      real(8), intent(in), optional :: atol
      logical                       :: flag

      flag = .true.

      if (present(rtol)) then
         flag = flag .and. (abs(actual - desired)/abs(desired)) <= rtol
      end if

      if (present(atol)) then
         flag = flag .and. (abs(actual - desired)) <= atol
      end if
   end function isclose_scalar

   !---------------------------------------------------------------------------
   ! isclose_vector
   !---------------------------------------------------------------------------
   function isclose_vector(actual, desired, rtol, atol) result( flag )
      real(8), dimension(:), intent(in)           :: actual
      real(8), dimension(:), intent(in)           :: desired
      real(8),               intent(in), optional :: rtol
      real(8),               intent(in), optional :: atol
      logical                                     :: flag

      flag = .true.

      if (present(rtol)) then
         flag = flag .and. maxval(abs(actual - desired)/abs(desired)) <= rtol
      end if

      if (present(atol)) then
         flag = flag .and. maxval(abs(actual - desired)) <= atol
      end if
   end function isclose_vector

end module isclose_module
