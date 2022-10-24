!==============================================================================
! MODULE: transformation_module
!
! Define the suite of available scalar transformations.
!
! Public procedures
! -----------------
! apply_transformation : apply to the input yielding the output.
!
! is_valid_transformation : check if the transformation is available.
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
module transformation_module
   implicit none
   private

   public :: apply_transformation
   public :: is_valid_transformation

contains

   !---------------------------------------------------------------------------
   ! FUNCTION: apply_transformation
   !
   ! Apply the designated scalar transformation.
   !
   ! Arguments
   ! ---------
   ! transformation : string
   !   transformation name
   !
   ! input : real scalar
   !   value to be transformed
   !
   ! Returns
   ! -------
   ! output : real scalar
   !   value of the transformed input
   !
   !---------------------------------------------------------------------------
   function apply_transformation(transformation, input) result(output)
      character(len=12), intent(in) :: transformation
      real(8),           intent(in) :: input
      real(8)                       :: output

      select case(transformation)
         case ("EXPONENTIAL")
            output = exp(input)

         case ("LINEAR")
            output = input

         case ("SOFTPLUS")
            output = log(exp(input) + 1.0_8)

         case ("TANH")
            output =  tanh(input)

         case default
            print *, transformation
            error stop "transformation type not implemented."

      end select
   end function apply_transformation

   !---------------------------------------------------------------------------
   ! FUNCTION: is_valid_transformation
   !
   ! Checks that the named transformation is available (valid).
   !
   ! Arguments
   ! ---------
   ! transformation : string
   !   transformation name
   !
   ! Returns
   ! -------
   ! is_valid : logical
   !   equals .true. if the transformation is available, otherwise .false.
   !
   !---------------------------------------------------------------------------
   function is_valid_transformation(transformation) result(is_valid)
      character(len=12), intent(in) :: transformation
      logical                       :: is_valid

      is_valid = .true.

      select case(transformation)
         case ("EXPONENTIAL")
         case ("LINEAR")
         case ("SOFTPLUS")
         case ("TANH")

         case default
            is_valid = .false.
      end select
   end function is_valid_transformation

end module transformation_module
