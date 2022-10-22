!==============================================================================
! MODULE: transformation_module
!
! Define the suite of available scalar transfomations.
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
! * 20 October 2022
!
!==============================================================================
module transformation_module
   implicit none
   private

   public :: apply_transformation
   public :: is_valid_transformation

contains

   !-----------------------------------
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
            stop "transformation type not implemented."

      end select
   end function apply_transformation

   !-----------------------------------
   function is_valid_transformation(transformation) result(is_valid)
      character(len=12), intent(in) :: transformation
      logical                       :: is_valid

      is_valid = .TRUE.

      select case(transformation)
         case ("EXPONENTIAL")
         case ("LINEAR")
         case ("SOFTPLUS")
         case ("TANH")

         case default
            is_valid = .FALSE.
      end select
   end function is_valid_transformation

end module transformation_module
