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

      x_sample = [ &
         4.0000e+00,  5.5000e+01, -2.6300e+01, -5.2000e+00, &
        -1.5800e+01,  4.7400e+01, -3.8900e+01,  9.4400e+01, &
        -3.8900e+01, -1.6700e+01,  1.4410e+02,  1.8800e+01, &
         4.7500e+01, -5.0000e+00,  1.8000e+01,  1.7800e+01, &
         2.6400e+01,  1.3978e+03,  5.0000e-01, -1.5000e+00, &
         5.0000e-01,  5.0000e-01 &
      ]

      desired = [ &
         -11.033484, -5.495101, 42.222622, 36.678978, 0.037328 &
      ]

      details = read_blueprint(filename)
      call initialize(test_model, details)
      actual = employ(test_model, x_sample)

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_model'
      end if
   end subroutine

end module test_model_module
