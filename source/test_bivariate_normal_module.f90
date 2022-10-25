!==============================================================================
! MODULE: test_bivariate_normal_module
!
! Test the Fortran-based bivariate normal distribution utility functions.
!
! Notes
! -----
! * The "desired" values were actual using MATLAB 2022b.
!
! * This code was tested using:
!   GNU Fortran (MinGW-W64 x86_64-ucrt-posix-seh, built by Brecht Sanders) 12.2.0
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
module test_bivariate_normal_module
   use bivariate_normal_module
   use isclose_module
   use iso_fortran_env, only : ERROR_UNIT

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: test_bivariate_normal

   !-----------------------------------
   ! Module parameters
   !-----------------------------------
   real(8), parameter :: ABSOLUTE_TOLERANCE = 1e-6
   real(8), parameter :: RELATIVE_TOLERANCE = 1e-4

   contains

   !-----------------------------------
   subroutine test_bivariate_normal()
      write(*,*) 'test_bivariate_normal'

      call test_bivariate_normal_pdf_elemental()
      call test_bivariate_normal_pdf_uv_array()

      call test_bivariate_normal_cdf_elemental()
      call test_bivariate_normal_cdf_uv_array()
   end subroutine test_bivariate_normal

   !-----------------------------------
   subroutine test_bivariate_normal_pdf_elemental()
      real(8), dimension(4) :: u, v, mu_u, mu_v, sigma_u, sigma_v, rho
      real(8), dimension( size(u) ) :: actual, desired

      u       = [8.0, 9.0, 2.0, 9.0]
      v       = [4.0, 9.0, 8.0, 9.0]
      mu_u    = [4.0, 4.0, 7.0, 8.0]
      mu_v    = [2.0, 5.0, 5.0, 6.0]
      sigma_u = [9.0, 4.0, 6.0, 3.0]
      sigma_v = [6.0, 2.0, 2.0, 5.0]
      rho     = [0.2, 0.3, 0.8, 0.3]

      actual  = pdf(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
      desired = [0.002641689, 0.002237701, 0.000023018, 0.009171173]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_bivariate_normal_pdf_elemental'
      end if
   end subroutine test_bivariate_normal_pdf_elemental

   !-----------------------------------
   subroutine test_bivariate_normal_pdf_uv_array()
      real(8), dimension(4) :: u, v
      real(8) :: mu_u, mu_v, sigma_u, sigma_v, rho
      real(8), dimension( size(u) ) :: actual, desired

      u       = [8.0, 9.0, 2.0, 9.0]
      v       = [4.0, 9.0, 8.0, 9.0]
      mu_u    = 4.0
      mu_v    = 2.0
      sigma_u = 9.0
      sigma_v = 6.0
      rho     = 0.2

      actual  = pdf(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
      desired = [0.002641689, 0.001442926, 0.001662721, 0.001442926]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_bivariate_normal_pdf_uv_array'
      end if
   end subroutine test_bivariate_normal_pdf_uv_array

   !-----------------------------------
   subroutine test_bivariate_normal_cdf_elemental()
      real(8), dimension(4) :: u, v, mu_u, mu_v, sigma_u, sigma_v, rho
      real(8), dimension( size(u) ) :: actual, desired

      u       = [8.0, 9.0, 2.0, 9.0]
      v       = [4.0, 9.0, 8.0, 9.0]
      mu_u    = [4.0, 4.0, 7.0, 8.0]
      mu_v    = [2.0, 5.0, 5.0, 6.0]
      sigma_u = [9.0, 4.0, 6.0, 3.0]
      sigma_v = [6.0, 2.0, 2.0, 5.0]
      rho     = [0.2, 0.3, 0.8, 0.3]

      actual  = cdf(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
      desired = [0.121805190516559, 0.89270174122853, 0.998958698546541, 0.175450574637629]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_bivariate_normal_cdf_elemental'
      end if
   end subroutine test_bivariate_normal_cdf_elemental

   !-----------------------------------
   subroutine test_bivariate_normal_cdf_uv_array()
      real(8), dimension(4) :: u, v
      real(8) :: mu_u, mu_v, sigma_u, sigma_v, rho
      real(8), dimension( size(u) ) :: actual, desired

      u       = [8.0, 9.0, 2.0, 9.0]
      v       = [4.0, 9.0, 8.0, 9.0]
      mu_u    = 4.0
      mu_v    = 2.0
      sigma_u = 9.0
      sigma_v = 6.0
      rho     = 0.2

      actual  = cdf(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
      desired = [0.121805190516559, 0.520318147269926, 0.447250185167959, 0.520318147269926]

      if (.not. isclose_vector(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_bivariate_normal_cdf_uv_array'
      end if
   end subroutine test_bivariate_normal_cdf_uv_array

end module test_bivariate_normal_module
