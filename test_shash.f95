!==============================================================================
! test_shash.f95                                                     (06.05.22)
!
! Test the FORTRAN-based sinh-arcsinh normal distribution utility functions.
!
! Notes
! -----
! * In all tests, the "truth" values for the test were computed using
!   TensorFlow [1].
!
! * This code was tested using:
!   GNU Fortran (x86_64-posix-seh-rev0, Built by MinGW-W64 project) 8.1.0
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
!==============================================================================

PROGRAM main
    CALL test_pdf_elemental()
    CALL test_pdf_x_array()

    CALL test_cdf()
    CALL test_quantile()
    CALL test_quantile_random()
    CALL test_median()
    CALL test_median_random()
END PROGRAM main


!------------------------------------------------------------------------------
SUBROUTINE test_pdf_elemental()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: x, mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: computed, truth

    x     = (/ 0.5, 1.0, 1.5, 2.0 /)
    mu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    sigma = (/ 0.5, 1.0, 1.5, 2.0 /)
    nu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    tau   = (/ 1.0, 1.0, 0.8, 1.5 /)

    computed = pdf(x, mu, sigma, nu, tau)
    truth    = (/ 0.48394145, 0.35682481, 0.21718473, 0.1577138 /)
    WRITE(*,*) 'test_pdf_elemental:        max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_pdf_elemental


!------------------------------------------------------------------------------
SUBROUTINE test_pdf_x_array()
    USE shash_module
    IMPLICIT NONE

    REAL(8) :: x(4), mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: computed, truth

    x     = (/ 1.5, 2.5, 3.5, 4.5 /)
    mu    = 1.5
    sigma = 2.0
    nu    = 1.5
    tau   = 1.5

    computed = pdf(x, mu, sigma, nu, tau)
    truth    = (/ 0.069731, 0.16752543, 0.11880247, 0.08460783 /)
    WRITE(*,*) 'test_pdf_x_array:          max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_pdf_x_array


!------------------------------------------------------------------------------
SUBROUTINE test_cdf()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: x, mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: computed, truth

    x     = (/ 0.5, 1.0, 1.5, 2.0 /)
    mu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    sigma = (/ 0.5, 1.0, 1.5, 2.0 /)
    nu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    tau   = (/ 1.0, 1.0, 0.8, 1.5 /)

    computed = cdf(x, mu, sigma, nu, tau)
    truth    = (/ 0.84134475, 0.4925046, 0.22218556, 0.07596765 /)
    WRITE(*,*) 'test_cdf:                  max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_cdf


!------------------------------------------------------------------------------
SUBROUTINE test_quantile()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: pr, mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(pr) ) :: computed, truth

    pr    = (/ 0.1, 0.3, 0.6, 0.8 /)
    mu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    sigma = (/ 0.5, 1.0, 1.5, 2.0 /)
    nu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    tau   = (/ 1.0, 1.0, 0.8, 1.5 /)

    computed = quantile(pr, mu, sigma, nu, tau)
    truth    = (/ -0.64077578, 0.49707086, 3.46847072, 15.376874 /)
    WRITE(*,*) 'test_quantile (x):         max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_quantile


!------------------------------------------------------------------------------
SUBROUTINE test_quantile_random()
    USE shash_module
    IMPLICIT NONE

    INTEGER, PARAMETER :: NTESTS = 1000
    REAL(8), DIMENSION(NTESTS) :: pr, mu, sigma, nu, tau
    REAL(8), DIMENSION(NTESTS) :: computed, computed_cdf

    CALL RANDOM_NUMBER(pr)
    CALL RANDOM_NUMBER(mu)
    CALL RANDOM_NUMBER(sigma)
    CALL RANDOM_NUMBER(nu)
    CALL RANDOM_NUMBER(tau)

    computed = quantile(pr, mu, sigma, nu, tau)
    computed_cdf = cdf(computed, mu, sigma, nu, tau)
    WRITE(*,*) 'test_quantile_random (p):  max_abs_error = ', MAXVAL(ABS(computed_cdf-pr))
END SUBROUTINE test_quantile_random


!------------------------------------------------------------------------------
SUBROUTINE test_median()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(mu) ) :: computed, truth

    mu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    sigma = (/ 0.5, 1.0, 1.5, 2.0 /)
    nu    = (/ 0.0, 0.5, 1.0, 1.5 /)
    tau   = (/ 1.0, 1.0, 0.8, 1.5 /)

    computed = median(mu, sigma, nu, tau)
    truth    = (/ 0.0, 1.0210953, 2.8640418, 5.8619227 /)
    WRITE(*,*) 'test_median (x):           max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_median


!------------------------------------------------------------------------------
SUBROUTINE test_median_random()
    USE shash_module
    IMPLICIT NONE

    INTEGER, PARAMETER :: NTESTS = 1000
    REAL(8), DIMENSION(NTESTS) :: pr, mu, sigma, nu, tau
    REAL(8), DIMENSION(NTESTS) :: computed, computed_cdf

    CALL RANDOM_NUMBER(pr)
    CALL RANDOM_NUMBER(mu)
    CALL RANDOM_NUMBER(sigma)
    CALL RANDOM_NUMBER(nu)
    CALL RANDOM_NUMBER(tau)

    computed = median(mu, sigma, nu, tau)
    computed_cdf = cdf(computed, mu, sigma, nu, tau)
    WRITE(*,*) 'test_median_random (p):    max_abs_error = ', MAXVAL(ABS(computed_cdf-0.5))
END SUBROUTINE test_median_random
