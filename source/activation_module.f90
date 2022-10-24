!==============================================================================
! MODULE: activation_module
!
! Define the suite of available layer activation functions.  This suite
! includes all activation functions in "Module:tf.keras.activations" of
! Tensorflow 2.10.
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
module activation_module
   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: apply_activation
   public :: is_valid_activation

contains

   !-----------------------------------
   function apply_activation(activation, input) result(output)
      character(len=12),              intent(in) :: activation
      real(8),          dimension(:), intent(in) :: input
      real(8),          dimension(size(input))   :: output

      ! SELU constants
      real(8), parameter :: ALPHA  = 1.6732632423543772848170429916717_8
      real(8), parameter :: LAMBDA = 1.0507009873554804934193349852946_8

      select case(trim(activation))
         case ("ELU")
            where (input > 0.0_8)
               output = input
            elsewhere
               output = exp(input) - 1.0_8
            end where

         case ("EXPONENTIAL")
            output = exp(input)

         case ("GELU")
            output = input * (1.0_8 + erf(input/sqrt(2.0_8))) / 2.0_8

         case ("HARD_SIGMOID")
            where (input < -2.5_8)
               output = 0.0_8
            elsewhere (input > 2.5_8)
               output = 1.0_8
            elsewhere
               output = 0.2_8 * input + 0.5_8
            end where

         case ("LINEAR")
            output = input

         case ("RELU")
            output = max(0.0_8, input)

         case ("SELU")
            where (input > 0.0_8)
               output = LAMBDA * input
            elsewhere
               output = LAMBDA * ALPHA * (exp(input) - 1.0_8)
            end where

         case ("SIGMOID")
            output = 1.0_8 / (1.0_8 + exp(-input))

         case ("SOFTMAX")
            output = exp(input)/sum(exp(input))

         case ("SOFTPLUS")
            output = log(exp(input) + 1.0_8)

         case ("SOFTSIGN")
            output = input / (abs(input) + 1.0_8)

         case ("SWISH")
            output =  input / (1.0_8 + exp(-input))

         case ("TANH")
            output =  tanh(input)

         case default
            stop "activation type not implemented."

      end select
   end function apply_activation

   !-----------------------------------
   function is_valid_activation(activation) result(is_valid)
      character(len=12), intent(in) :: activation
      logical                       :: is_valid

      is_valid = .TRUE.

      select case(activation)
         case ("ELU")
         case ("EXPONENTIAL")
         case ("GELU")
         case ("HARD_SIGMOID")
         case ("LINEAR")
         case ("RELU")
         case ("SELU")
         case ("SIGMOID")
         case ("SOFTMAX")
         case ("SOFTPLUS")
         case ("SOFTSIGN")
         case ("SWISH")
         case ("TANH")

         case default
            is_valid = .FALSE.
      end select
   end function is_valid_activation

end module activation_module
