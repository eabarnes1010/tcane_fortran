!==============================================================================
! MODULE: test_activation_module
!
! Test the available activation functions.
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
module test_activation_module
   use activation_module
   use isclose_module
   use iso_fortran_env, only : ERROR_UNIT

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: test_activation

   !-----------------------------------
   ! Module parameters
   !-----------------------------------
   real(8), parameter :: ABSOLUTE_TOLERANCE = 1e-5

   contains

   !-----------------------------------
   ! SUBROUTINE: test_activation
   !-----------------------------------
   subroutine test_activation()
      write(*,*) 'test_activation'

      call test_apply_activation
   end subroutine

   !-----------------------------------
   ! SUBROUTINE: test_apply_activation
   !-----------------------------------
   subroutine test_apply_activation()
      real(8), dimension(5) :: input  = [-3.0, -1.0, 1.0, 2.0, 3.0]
      real(8), dimension(5) :: actual, desired

      actual  = apply_activation('ELU         ', input)
      desired = [-0.950212931632136, -0.632120558828558, 1.0, 2.0, 3.0]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_ELU'
      end if

      actual  = apply_activation('EXPONENTIAL ', input)
      desired = [0.049787068367864, 0.367879441171442, 2.718281828459045, 7.389056098930650, 20.085536923187668]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_EXPONENTIAL'
      end if

      actual  = apply_activation('GELU        ', input)
      desired = [-0.004049694094890, -0.158655253931457, 0.841344746068543, 1.954499736103642, 2.995950305905110]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_GELU'
      end if

      actual  = apply_activation('HARD_SIGMOID', input)
      desired = [0.0, 0.3, 0.7, 0.9, 1.0]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_HARD_SIGMOID'
      end if

      actual  = apply_activation('LINEAR      ', input)
      desired = [-3.0, -1.0, 1.0, 2.0, 3.0]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_LINEAR'
      end if

      actual  = apply_activation('RELU        ', input)
      desired = [0.0, 0.0, 1.0, 2.0, 3.0]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_RELU'
      end if

      actual  = apply_activation('SELU        ', input)
      desired = [-1.670568728767112, -1.111330737812563, 1.050700987355480, 2.101401974710961, 3.152102962066441]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_SELU'
      end if

      actual  = apply_activation('SIGMOID     ', input)
      desired = [0.047425873177567, 0.268941421369995, 0.731058578630005, 0.880797077977882, 0.952574126822433]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_SIGMOID'
      end if

      actual  = apply_activation('SOFTMAX     ', input)
      desired = [0.001626468077848, 0.012018063870336, 0.088802148138444, 0.241389265612860, 0.656164054300512]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_SOFTMAX'
      end if

      actual  = apply_activation('SOFTPLUS    ', input)
      desired = [0.048587351573742, 0.313261687518223, 1.313261687518223, 2.126928011042973, 3.048587351573742]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_SOFTPLUS'
      end if

      actual  = apply_activation('SOFTSIGN    ', input)
      desired = [-0.75, -0.50, 0.50, 0.666666666666667, 0.75]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_SOFTSIGN'
      end if

      actual  = apply_activation('SWISH       ', input)
      desired = [-0.142277619532700, -0.268941421369995, 0.731058578630005, 1.761594155955765, 2.857722380467300]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_SWISH'
      end if

      actual  = apply_activation('TANH        ', input)
      desired = [-0.995054753686730, -0.761594155955765, 0.761594155955765, 0.964027580075817, 0.995054753686730]
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_activation_TANH'
      end if
   end subroutine test_apply_activation

end module test_activation_module
