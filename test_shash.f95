!==============================================================================
! test_shash.f95                                                     (14.10.22)
!
! Test the FORTRAN-based sinh-arcsinh normal (SHASH) distribution utility
! functions.
!
! Notes
! -----
! * The pdf and cdf "truth" values were computed using TensorFlow [1]. The
!   other "truth" values were computed using the python implementation.
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
! * 14 October 2022
!
!==============================================================================
SUBROUTINE test_shash()
    WRITE(*,*) '=========='
    WRITE(*,*) 'test_shash'
    WRITE(*,*) '=========='

    CALL test_shash_pdf_scalar()
    CALL test_shash_pdf_elemental()
    CALL test_shash_pdf_x_array()

    CALL test_shash_cdf_scalar()
    CALL test_shash_cdf_elemental()
    CALL test_shash_cdf_x_array()

    CALL test_shash_quantile()
    CALL test_shash_quantile_random()

    CALL test_shash_median()
    CALL test_shash_median_random()

    CALL test_shash_mean()
    CALL test_shash_mode()
    CALL test_shash_variance()
    CALL test_shash_stddev()
    CALL test_shash_skew()

    CALL test_shash_summary_statistics()
END SUBROUTINE test_shash


!------------------------------------------------------------------------------
SUBROUTINE test_shash_pdf_scalar()
    USE shash_module
    IMPLICIT NONE

    REAL(8) :: x, mu, sigma, nu, tau
    REAL(8) :: computed, truth

    x     = 4.5
    mu    = 1.5
    sigma = 2.0
    nu    = 1.5
    tau   = 1.5

    computed = pdf(x, mu, sigma, nu, tau)
    truth    = 0.08460783
    WRITE(*,*) 'test_shash_pdf_scalar:           max_abs_error = ', ABS(computed-truth)
END SUBROUTINE test_shash_pdf_scalar


!------------------------------------------------------------------------------
SUBROUTINE test_shash_pdf_elemental()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: x, mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: computed, truth

    x     = [0.5, 1.0, 1.5, 2.0]
    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = pdf(x, mu, sigma, nu, tau)
    truth    = [0.48394145, 0.35682481, 0.21718473, 0.1577138]
    WRITE(*,*) 'test_shash_pdf_elemental:        max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_pdf_elemental


!------------------------------------------------------------------------------
SUBROUTINE test_shash_pdf_x_array()
    USE shash_module
    IMPLICIT NONE

    REAL(8) :: x(4), mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: computed, truth

    x     = [1.5, 2.5, 3.5, 4.5]
    mu    = 1.5
    sigma = 2.0
    nu    = 1.5
    tau   = 1.5

    computed = pdf(x, mu, sigma, nu, tau)
    truth    = [0.069731, 0.16752543, 0.11880247, 0.08460783]
    WRITE(*,*) 'test_shash_pdf_x_array:          max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_pdf_x_array


!------------------------------------------------------------------------------
SUBROUTINE test_shash_cdf_scalar()
    USE shash_module
    IMPLICIT NONE

    REAL(8) :: x, mu, sigma, nu, tau
    REAL(8) :: computed, truth

    x     = 1.5
    mu    = 1.0
    sigma = 1.5
    nu    = 1.0
    tau   = 0.8

    computed = cdf(x, mu, sigma, nu, tau)
    truth    = 0.22218556
    WRITE(*,*) 'test_shash_cdf_scalar:           max_abs_error = ', ABS(computed-truth)
END SUBROUTINE test_shash_cdf_scalar


!------------------------------------------------------------------------------
SUBROUTINE test_shash_cdf_elemental()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: x, mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: computed, truth

    x     = [0.5, 1.0, 1.5, 2.0]
    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = cdf(x, mu, sigma, nu, tau)
    truth    = [0.84134475, 0.4925046, 0.22218556, 0.07596765]
    WRITE(*,*) 'test_shash_cdf_elemental:        max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_cdf_elemental


!------------------------------------------------------------------------------
SUBROUTINE test_shash_cdf_x_array()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: x, mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: computed, truth

    x     = [1.5, 2.5, 3.5, 4.5]
    mu    = 1.5
    sigma = 2.0
    nu    = 1.5
    tau   = 1.5

    computed = cdf(x, mu, sigma, nu, tau)
    truth    = [0.01661558, 0.1599833, 0.3035546, 0.4036644]
    WRITE(*,*) 'test_shash_cdf_x_array:          max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_cdf_x_array


!------------------------------------------------------------------------------
SUBROUTINE test_shash_quantile()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: pr, mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(pr) ) :: computed, truth

    pr    = [0.1, 0.3, 0.6, 0.8]
    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = quantile(pr, mu, sigma, nu, tau)
    truth    = [-0.64077578, 0.49707086, 3.46847072, 15.376874]
    WRITE(*,*) 'test_shash_quantile (x):         max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_quantile


!------------------------------------------------------------------------------
SUBROUTINE test_shash_quantile_random()
    USE shash_module
    IMPLICIT NONE

    INTEGER, PARAMETER :: NTESTS = 1000
    REAL(8), DIMENSION(NTESTS) :: pr, mu, sigma, nu, tau
    REAL(8), DIMENSION(NTESTS) :: computed, computed_cdf

    CALL RANDOM_NUMBER(pr)
    CALL RANDOM_NUMBER(mu)
    mu = 2.0*(mu - 0.5)
    CALL RANDOM_NUMBER(sigma)
    sigma = 2.0*sigma
    CALL RANDOM_NUMBER(nu)
    nu = 2.0*(nu - 0.5)
    CALL RANDOM_NUMBER(tau)
    tau = 2.0*tau

    computed = quantile(pr, mu, sigma, nu, tau)
    computed_cdf = cdf(computed, mu, sigma, nu, tau)
    WRITE(*,*) 'test_shash_quantile_random (p):  max_abs_error = ', MAXVAL(ABS(computed_cdf-pr))
END SUBROUTINE test_shash_quantile_random


!------------------------------------------------------------------------------
SUBROUTINE test_shash_median()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(mu) ) :: computed, truth

    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = median(mu, sigma, nu, tau)
    truth    = [0.0, 1.0210953, 2.8640418, 5.8619227]
    WRITE(*,*) 'test_shash_median (x):           max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_median


!------------------------------------------------------------------------------
SUBROUTINE test_shash_median_random()
    USE shash_module
    IMPLICIT NONE

    INTEGER, PARAMETER :: NTESTS = 1000
    REAL(8), DIMENSION(NTESTS) :: mu, sigma, nu, tau
    REAL(8), DIMENSION(NTESTS) :: computed, computed_cdf

    CALL RANDOM_NUMBER(mu)
    mu = 2.0*(mu - 0.5)
    CALL RANDOM_NUMBER(sigma)
    sigma = 2.0*sigma
    CALL RANDOM_NUMBER(nu)
    nu = 2.0*(nu - 0.5)
    CALL RANDOM_NUMBER(tau)
    tau = 2.0*tau

    computed = median(mu, sigma, nu, tau)
    computed_cdf = cdf(computed, mu, sigma, nu, tau)
    WRITE(*,*) 'test_shash_median_random (p):    max_abs_error = ', MAXVAL(ABS(computed_cdf-0.5))
END SUBROUTINE test_shash_median_random


!------------------------------------------------------------------------------
SUBROUTINE test_shash_mean()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(mu) ) :: computed, truth

    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = mean(mu, sigma, nu, tau)
    truth    = [0.0, 1.2058396, 3.2700596, 9.874768]
    WRITE(*,*) 'test_shash_mean:                 max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_mean


!------------------------------------------------------------------------------
SUBROUTINE test_shash_mode()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(mu) ) :: computed, truth

    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = mode(mu, sigma, nu, tau)
    truth    = [0.0, 0.5764171, 1.6523985, 2.304149]
    WRITE(*,*) 'test_shash_mode:                 max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_mode


!------------------------------------------------------------------------------
SUBROUTINE test_shash_variance()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(mu) ) :: computed, truth

    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = variance(mu, sigma, nu, tau)
    truth    = [0.25, 1.3164116, 4.4789896, 105.18572]
    WRITE(*,*) 'test_shash_variance:             max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_variance


!------------------------------------------------------------------------------
SUBROUTINE test_shash_stddev()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(mu) ) :: computed, truth

    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = stddev(mu, sigma, nu, tau)
    truth    = [0.5, 1.1473497, 2.1163623, 10.256009]
    WRITE(*,*) 'test_shash_stddev:               max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_stddev


!------------------------------------------------------------------------------
SUBROUTINE test_shash_skew()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(mu) ) :: computed, truth

    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    computed = skew(mu, sigma, nu, tau)
    truth    = [0.0, 0.7544219, 0.7953873, 2.29442]
    WRITE(*,*) 'test_shash_skew:                 max_abs_error = ', MAXVAL(ABS(computed-truth))
END SUBROUTINE test_shash_skew


!------------------------------------------------------------------------------
SUBROUTINE test_shash_summary_statistics()
    USE shash_module
    IMPLICIT NONE

    REAL(8), DIMENSION(4) :: mu, sigma, nu, tau
    TYPE (SUMMARY_STATISTICS), DIMENSION(4) :: computed, truth

    mu    = [0.0, 0.5, 1.0, 1.5]
    sigma = [0.5, 1.0, 1.5, 2.0]
    nu    = [0.0, 0.5, 1.0, 1.5]
    tau   = [1.0, 1.0, 0.8, 1.5]

    WRITE(*,*) 'test_shash_summary_statistics'
END SUBROUTINE test_shash_summary_statistics
