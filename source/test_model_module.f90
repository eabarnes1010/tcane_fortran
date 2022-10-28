!==============================================================================
! MODULE: test_model_module
!
! Test the Fortran-based evaluation system for the Shash and Bivariate Normal
! Artificial Neural Network functions.
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
module test_model_module
   use blueprint_module
   use isclose_module
   use iso_fortran_env, only : ERROR_UNIT
   use model_module

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: test_model

   !-----------------------------------
   ! Module parameters
   !-----------------------------------
   real(8), parameter :: ABSOLUTE_TOLERANCE = 1e-5

   contains

   !-----------------------------------
   ! SUBROUTINE: test_model
   !-----------------------------------
   subroutine test_model()
      write(*,*) 'test_model'

      call simple_test()
   end subroutine

   !-----------------------------------
   ! SUBROUTINE: simple_test
   !-----------------------------------
   subroutine simple_test()
      real(8), dimension(22) :: x_sample
      real(8), dimension(5)  :: actual, desired
      type(Model) :: test_model
      type(Blueprint) :: details
      character(len=*), parameter :: filename = ".\data\test_blueprint.json"

      details = read_blueprint(filename)
      call initialize(test_model, details)

      x_sample = [ &
           4.00,  50.00,  74.70, -10.70, -32.00, -32.00,  -8.30, -30.60, &
          47.20,  -8.30, 111.40,  16.30,  63.20,  10.00,  16.00,  15.70, &
          29.10, 654.20,   9.80,   3.80,  -2.20, -11.20 &
      ]

      desired = [ &
         18.785368, -5.833050, 38.011936, 38.096218, -0.021327 &
      ]

      actual = employ(test_model, x_sample)
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_model'
      end if

      x_sample = [ &
           4.00,  125.00,  13.40, -29.50,  -8.00,  24.10,  16.70, -16.70, &
         -16.70,   16.70, 127.60,  15.40, 102.00,  10.00,  14.70,   6.90, &
          26.70, 1815.20,   0.00,  -8.00,   4.00,   4.00 &
      ]

      desired = [ &
         2.902810, -13.442649, 32.559185, 30.183308, 0.029761 &
      ]

      actual = employ(test_model, x_sample)
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_model'
      end if
   end subroutine

end module test_model_module
