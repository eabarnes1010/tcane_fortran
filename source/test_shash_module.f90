!==============================================================================
! MODULE: test_shash_module
!
! Test the Fortran-based sinh-arcsinh normal (SHASH) distribution utility
! functions.
!
! Notes
! -----
! * The pdf and cdf "desired" values were actual using TensorFlow [1]. The
!   other "desired" values were actual using the python implementation.
!
! * This code was tested using:
!   GNU Fortran (MinGW-W64 x86_64-ucrt-posix-seh, built by Brecht Sanders) 12.2.0
!
! References
! ----------
! [1] tfp.distributions.SinhArcsinh. url: https://www.tensorflow.org/
! probability/api_docs/python/tfp/distributions/SinhArcsinh.
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
module test_shash_module
   use isclose_module
   use iso_fortran_env, only : ERROR_UNIT
   use shash_module

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: test_shash

   !-----------------------------------
   ! Module parameters
   !-----------------------------------
   real(8), parameter :: ABSOLUTE_TOLERANCE = 1e-6
   real(8), parameter :: RELATIVE_TOLERANCE = 1e-4

   contains

   !-----------------------------------
   subroutine test_shash()
      write(*,*) 'test_shash'

      call test_shash_pdf_elemental()
      call test_shash_pdf_x_array()

      call test_shash_cdf_elemental()
      call test_shash_cdf_x_array()

      call test_shash_quantile()
      call test_shash_quantile_random()

      call test_shash_median()
      call test_shash_median_random()

      call test_shash_mean()
      call test_shash_mode()
      call test_shash_variance()
      call test_shash_stddev()
      call test_shash_skew()
   end subroutine test_shash

   !-----------------------------------
   subroutine test_shash_pdf_elemental()
      real(8), dimension(4) :: x, loc, sigma, skewness, tailweight
      real(8), dimension( size(x) ) :: actual, desired

      x          = [0.5, 1.0, 1.5, 2.0]
      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = pdf(x, loc, sigma, skewness, tailweight)
      desired = [0.48394145, 0.35682481, 0.21718473, 0.1577138]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_pdf_elemental'
      end if
   end subroutine test_shash_pdf_elemental

   !-----------------------------------
   subroutine test_shash_pdf_x_array()
      real(8) :: x(4), loc, sigma, skewness, tailweight
      real(8), dimension( size(x) ) :: actual, desired

      x          = [1.5, 2.5, 3.5, 4.5]
      loc        = 1.5
      sigma      = 2.0
      skewness   = 1.5
      tailweight = 1.5

      actual  = pdf(x, loc, sigma, skewness, tailweight)
      desired = [0.069731, 0.16752543, 0.11880247, 0.08460783]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_pdf_x_array'
      end if
   end subroutine test_shash_pdf_x_array

   !-----------------------------------
   subroutine test_shash_cdf_elemental()
      real(8), dimension(4) :: x, loc, sigma, skewness, tailweight
      real(8), dimension( size(x) ) :: actual, desired

      x          = [0.5, 1.0, 1.5, 2.0]
      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = cdf(x, loc, sigma, skewness, tailweight)
      desired = [0.84134475, 0.4925046, 0.22218556, 0.07596765]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_cdf_elemental'
      end if
   end subroutine test_shash_cdf_elemental

   !-----------------------------------
   subroutine test_shash_cdf_x_array()
      real(8), dimension(4) :: x, loc, sigma, skewness, tailweight
      real(8), dimension( size(x) ) :: actual, desired

      x          = [1.5, 2.5, 3.5, 4.5]
      loc        = 1.5
      sigma      = 2.0
      skewness   = 1.5
      tailweight = 1.5

      actual  = cdf(x, loc, sigma, skewness, tailweight)
      desired = [0.01661558, 0.1599833, 0.3035546, 0.4036644]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_cdf_x_array'
      end if
   end subroutine test_shash_cdf_x_array

   !-----------------------------------
   subroutine test_shash_quantile()
      real(8), dimension(4) :: pr, loc, sigma, skewness, tailweight
      real(8), dimension( size(pr) ) :: actual, desired

      pr         = [0.1, 0.3, 0.6, 0.8]
      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = quantile(pr, loc, sigma, skewness, tailweight)
      desired = [-0.64077578, 0.49707086, 3.46847072, 15.376874]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_quantile'
      end if
   end subroutine test_shash_quantile

   !-----------------------------------
   subroutine test_shash_quantile_random()
      integer, parameter :: NTESTS = 1000
      real(8), dimension(NTESTS) :: pr, loc, sigma, skewness, tailweight
      real(8), dimension(NTESTS) :: actual, q

      call random_number(pr)
      call random_number(loc)
      loc = 2.0*(loc - 0.5)
      call random_number(sigma)
      sigma = 2.0*sigma
      call random_number(skewness)
      skewness = 2.0*(skewness - 0.5)
      call random_number(tailweight)
      tailweight = 2.0*tailweight

      q      = quantile(pr, loc, sigma, skewness, tailweight)
      actual = cdf(q, loc, sigma, skewness, tailweight)

      if (.not. isclose(actual, pr, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_quantile_random'
      end if
   end subroutine test_shash_quantile_random

   !-----------------------------------
   subroutine test_shash_median()
      real(8), dimension(4) :: loc, sigma, skewness, tailweight
      real(8), dimension( size(loc) ) :: actual, desired

      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = median(loc, sigma, skewness, tailweight)
      desired = [0.0, 1.0210953, 2.8640418, 5.8619227]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_median'
      end if
   end subroutine test_shash_median

   !-----------------------------------
   subroutine test_shash_median_random()
      integer, parameter :: NTESTS = 1000
      real(8), dimension(NTESTS) :: loc, sigma, skewness, tailweight
      real(8), dimension(NTESTS) :: actual, desired, med

      call random_number(loc)
      loc = 2.0*(loc - 0.5)
      call random_number(sigma)
      sigma = 2.0*sigma
      call random_number(skewness)
      skewness = 2.0*(skewness - 0.5)
      call random_number(tailweight)
      tailweight = 2.0*tailweight

      med     = median(loc, sigma, skewness, tailweight)
      actual  = cdf(med, loc, sigma, skewness, tailweight)
      desired = 0.5_8

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_median_random'
      end if
   end subroutine test_shash_median_random

   !-----------------------------------
   subroutine test_shash_mean()
      real(8), dimension(4) :: loc, sigma, skewness, tailweight
      real(8), dimension( size(loc) ) :: actual, desired

      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = mean(loc, sigma, skewness, tailweight)
      desired = [0.0, 1.2058396, 3.2700596, 9.874768]

      if (.not. isclose(actual, desired, atol=RELATIVE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_mean'
      end if
   end subroutine test_shash_mean

   !-----------------------------------
   subroutine test_shash_mode()
      real(8), dimension(4) :: loc, sigma, skewness, tailweight
      real(8), dimension( size(loc) ) :: actual, desired

      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = mode(loc, sigma, skewness, tailweight)
      desired = [0.0, 0.5764171, 1.6523985, 2.304149]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_mode'
      end if
   end subroutine test_shash_mode

   !-----------------------------------
   subroutine test_shash_variance()
      real(8), dimension(4) :: loc, sigma, skewness, tailweight
      real(8), dimension( size(loc) ) :: actual, desired

      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = variance(loc, sigma, skewness, tailweight)
      desired = [0.25, 1.3164116, 4.4789896, 105.18572]

      if (.not. isclose(actual, desired, atol=RELATIVE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_variance'
      end if
   end subroutine test_shash_variance

   !-----------------------------------
   subroutine test_shash_stddev()
      real(8), dimension(4) :: loc, sigma, skewness, tailweight
      real(8), dimension( size(loc) ) :: actual, desired

      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = stddev(loc, sigma, skewness, tailweight)
      desired = [0.5, 1.1473497, 2.1163623, 10.256009]

      if (.not. isclose(actual, desired, atol=RELATIVE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_stddev'
      end if
   end subroutine test_shash_stddev

   !-----------------------------------
   subroutine test_shash_skew()
      real(8), dimension(4) :: loc, sigma, skewness, tailweight
      real(8), dimension( size(loc) ) :: actual, desired

      loc        = [0.0, 0.5, 1.0, 1.5]
      sigma      = [0.5, 1.0, 1.5, 2.0]
      skewness   = [0.0, 0.5, 1.0, 1.5]
      tailweight = [1.0, 1.0, 0.8, 1.5]

      actual  = skew(loc, sigma, skewness, tailweight)
      desired = [0.0, 0.7544219, 0.7953873, 2.29442]

      if (.not. isclose(actual, desired, atol=ABSOLUTE_TOLERANCE)) then
         write(ERROR_UNIT,*) 'TEST FAILED: test_shash_skew'
      end if
   end subroutine test_shash_skew

end module test_shash_module
