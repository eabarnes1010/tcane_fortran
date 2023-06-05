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
! * 03 November 2022
!
!==============================================================================
module test_model_module
   use blueprint_module
   use model_module

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: test_model

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
      type(Model) :: test_model
      type(Blueprint) :: details

      details = read_blueprint("./data/intensity_shash3_test_blueprint.json")
      call initialize(test_model, details)
   end subroutine test_shash3_case

   !-----------------------------------
   ! SUBROUTINE: test_bivariate_normal_case
   !-----------------------------------
   subroutine test_bivariate_normal_case()
      type(Model) :: test_model
      type(Blueprint) :: details

      details = read_blueprint("./data/track_bivariate_normal_test_blueprint.json")
      call initialize(test_model, details)
   end subroutine test_bivariate_normal_case

end module test_model_module
