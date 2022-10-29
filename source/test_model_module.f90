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

      call test_shash3_case()
      call test_bivariate_normal_case()
   end subroutine

   !-----------------------------------
   ! SUBROUTINE: test_shash3_case
   !-----------------------------------
   subroutine test_shash3_case()
      real(8), dimension(16) :: x_sample
      real(8), dimension(4)  :: actual, desired
      type(Model) :: test_model
      type(Blueprint) :: details

      details = read_blueprint(".\data\intensity_shash3_test_blueprint.json")
      call initialize(test_model, details)

      x_sample = [ &
         30.0, 43.5,   4.0,   1.5,  4.5, -5.5, 1.5, -0.5, 0.0, 11.1, &
         28.6,  7.1, 812.4, -53.5, 97.0, 75.0 &
      ]

      desired = [ &
         0.23195606470108032, 9.098888397216797, -0.18818049132823944, 1.0 &
      ]

      actual = employ(test_model, x_sample)
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash3_case'
      end if
   end subroutine test_shash3_case

   !-----------------------------------
   ! SUBROUTINE: test_bivariate_normal_case
   !-----------------------------------
   subroutine test_bivariate_normal_case()
      real(8), dimension(22) :: x_sample
      real(8), dimension(5)  :: actual, desired
      type(Model) :: test_model
      type(Blueprint) :: details

      details = read_blueprint(".\data\track_bivariate_normal_test_blueprint.json")
      call initialize(test_model, details)

      x_sample = [ &
           4.0, 55.0, 28.7, -13.0, 28.5, -44.2, -52.8, -75.0, 125.0, 2.8, &
         111.1, 20.5, 37.0,   5.0, 25.0,   9.4,  30.3,  11.6,   3.0, 7.0, &
          -5.0, -5.0 &
      ]

      desired = [ &
         20.35586929321289, 2.1636693477630615, 46.34773635864258, &
         45.1985969543457,  0.008352618664503098 &
      ]

      actual = employ(test_model, x_sample)
      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_bivariate_normal_case'
      end if
   end subroutine test_bivariate_normal_case

end module test_model_module
