!==============================================================================
! shash.f95                                                          (06.04.22)
!
! FORTRAN-based sinh-arcsinh normal (SHASH) distribution utility functions.
!
! Public functions
! ----------------
! pdf(x, mu, sigma, nu, tau)
!   Compute the SHASH probability density function (pdf).
!
! cdf(x, mu, sigma, nu, tau)
!   Compute the SHASH cumulative distribution function (cdf).
!
! quantile(pr, mu, sigma, nu, tau)
!   Compute the SHASH inverse cumulative distribution function: that is,
!   find x such that cdf(x) = pr. This is an approximate function.
!
! median(mu, sigma, nu, tau)
!   Compute the distribution median.
!
! Notes
! -----
! * The sinh-arcsinh normal distribution was defined in [1]. A more accessible,
! though less comprehensive, presentation is given in [2].
!
! * These functions are all written as ELEMENTAL.  As such, they may be called
! with scalar arguments, or with array arguments -- as long as all of the
! arrays are commensurate.
!
! * This code was tested using:
!   GNU Fortran (x86_64-posix-seh-rev0, Built by MinGW-W64 project) 8.1.0
!
! References
! ----------
! [1] Jones, M. C. & Pewsey, A., Sinh-arcsinh distributions,
! Biometrika, Oxford University Press, 2009, 96, 761-780.
! DOI: 10.1093/biomet/asp053.
!
! [2] Jones, C. & Pewsey, A., The sinh-arcsinh normal distribution,
! Significance, Wiley, 2019, 16, 6-7.
! DOI: 10.1111/j.1740-9713.2019.01245.x.
! https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1740-9713.2019.01245.x
!
! [3] Abramowitz, M. & Stegun, I. A., Handbook of Mathematical Functions:
! with Formulas, Graphs, and Mathematical Tables, Dover Publications, 1965.
! ISBN: 978-0486612720.
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
MODULE shash_module
IMPLICIT NONE

PUBLIC  :: pdf, cdf, quantile, median

PRIVATE :: convert_tf_to_jp, rational_approximation, gaussian_cdf_inv

!---------------------------------------------------------------------------
! Mathematical constants
!---------------------------------------------------------------------------
REAL(8), PARAMETER, PRIVATE :: ONE_OVER_SQRT_TWO    = 0.7071067811865475244008444
REAL(8), PARAMETER, PRIVATE :: ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399461

CONTAINS

!------------------------------------------------------------------------------
! FUNCTION pdf(x, mu, sigma, nu, tau)
!
! Compute the SHASH probability density function (pdf).
!
! Parameters
! ----------
! x : real scalar
!   The value at which to compute the SHASH probability density function.
!
! mu : real scalar
!   The location parameter.
!   TensorFlow calls this "loc".
!
! sigma : real scalar
!   The scale parameter. Must be strictly positive.
!   TensorFlow calls this "scale".
!
! nu : real scalar
!   The skewness parameter.
!   TensorFlow calls this "skewness".
!
! tau : real scalar
!   The tail-weight parameter. Must be strictly positive.
!   TensorFlow calls this "tailweight".
!
! Returns
! -------
! pdf : real scalar
!   The computed shash(mu, sigma, nu, tau) probability density function
!   evaluated at x.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION pdf(x, mu, sigma, nu, tau)
    REAL(8) :: pdf
    REAL(8), INTENT(IN) :: x, mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: y, z

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) Equation (2) on page 762.
    y = (x - xi) / eta
    z = SINH(delta * ASINH(y) - eps)
    pdf = ONE_OVER_SQRT_TWO_PI * (delta / eta) * SQRT((1.0 + z*z) / (1.0 + y*y)) * EXP(-z*z / 2.0)
END FUNCTION pdf


!------------------------------------------------------------------------------
! FUNCTION cdf(x, mu, sigma, nu, tau)
!
! Compute the SHASH cumulative distribution function (cdf).
!
! Parameters
! ----------
! x : real scalar
!   The value at which to compute the SHASH cumulative distribution function.
!
! mu : real scalar
!   The location parameter.
!   TensorFlow calls this "loc".
!
! sigma : real scalar
!   The scale parameter. Must be strictly positive.
!   TensorFlow calls this "scale".
!
! nu : real scalar
!   The skewness parameter.
!   TensorFlow calls this "skewness".
!
! tau : real scalar
!   The tail-weight parameter. Must be strictly positive.
!   TensorFlow calls this "tailweight".
!
! Returns
! -------
! cdf : real scalar
!   The computed shash(mu, sigma, nu, tau) cumulative distribution function
!   evaluated at x.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION cdf(x, mu, sigma, nu, tau)
    REAL(8) :: cdf
    REAL(8), INTENT(IN) :: x, mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: y, z

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762.
    y = (x - xi) / eta
    z = SINH(delta * ASINH(y) - eps)
    cdf = 0.5 * (1.0 + ERF(ONE_OVER_SQRT_TWO * z))
END FUNCTION cdf


!------------------------------------------------------------------------------
! FUNCTION quantile(pr, mu, sigma, nu, tau)
!
! Compute the SHASH inverse cumulative distribution function: that is,
! find x such that cdf(x) = pr.
!
! Parameters
! ----------
! pr : real scalar in the range 0 < p < 1.
!   The probability value at which to compute the SHASH inverse cumulative
!   distribution function.
!
! mu : real scalar
!   The location parameter.
!   TensorFlow calls this "loc".
!
! sigma : real scalar
!   The scale parameter. Must be strictly positive.
!   TensorFlow calls this "scale".
!
! nu : real scalar
!   The skewness parameter.
!   TensorFlow calls this "skewness".
!
! tau : real scalar
!   The tail-weight parameter. Must be strictly positive.
!   TensorFlow calls this "tailweight".
!
! Returns
! -------
! quantile : real scalar
!   The computed shash(mu, sigma, nu, tau) inverse cumulative distribution
!   function evaluated at probability pr.
!
! Notes
! -----
! * This function uses Newton's method, with a pretty good starting guess.
!   The result is accurate but a bit slow.
!
! * If pr is out of the range 0 < pr < 1, this function fails.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION quantile(pr, mu, sigma, nu, tau)
    REAL(8) :: quantile
    REAL(8), INTENT(IN) :: pr, mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: z, x

    INTEGER :: iteration
    INTEGER, PARAMETER :: MAX_ITERATION = 5

    REAL(8) :: difference
    REAL(8), PARAMETER :: MAX_ABS_DIFFERENCE = 1.0e-6

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762, and an approximate
    ! Gaussian cdf inverse, to get an initial approximation.
    z = gaussian_cdf_inv(pr)
    x = xi + eta * SINH((ASINH(z) + eps) / delta)

    ! Improve the approximation using Newton's method.
    DO iteration = 1, MAX_ITERATION
        difference = cdf(x, mu, sigma, nu, tau) - pr
        x = x - difference/pdf(x, mu, sigma, nu, tau)
        IF (ABS(difference) <= MAX_ABS_DIFFERENCE) EXIT
    END DO

    quantile = x
END FUNCTION quantile


!------------------------------------------------------------------------------
! FUNCTION median(mu, sigma, nu, tau)
!
! Compute the distribution median.
!
! Parameters
! ----------
! mu : real scalar
!   The location parameter.
!   TensorFlow calls this "loc".
!
! sigma : real scalar
!   The scale parameter. Must be strictly positive.
!   TensorFlow calls this "scale".
!
! nu : real scalar
!   The skewness parameter.
!   TensorFlow calls this "skewness".
!
! tau : real scalar
!   The tail-weight parameter. Must be strictly positive.
!   TensorFlow calls this "tailweight".
!
! Returns
! -------
! median : real scalar
!   The computed median of shash(mu, sigma, nu, tau).
!
! Notes
! -----
! * Unlike the quantile function, this is an exact computation.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION median(mu, sigma, nu, tau)
    REAL(8) :: median
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762.
    median = xi + eta * SINH(eps / delta)
END FUNCTION median


!==============================================================================
! PRIVATE HELPER FUNCTIONS
!==============================================================================

!------------------------------------------------------------------------------
! SUBROUTINE convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)
!
! Convert from TensorFlow formulation and notation to the Jones and Pewsey
! formulation and notation.
!
! Parameters
! ----------
! mu : real scalar
!   The location parameter.
!   TensorFlow calls this "loc".
!
! sigma : real scalar
!   The scale parameter. Must be strictly positive.
!   TensorFlow calls this "scale".
!
! nu : real scalar
!   The skewness parameter.
!   TensorFlow calls this "skewness".
!
! tau : real scalar
!   The tail-weight parameter. Must be strictly positive.
!   TensorFlow calls this "tailweight".
!
! Notes
! -----
! * This conversion of the parameters is necessary to use the detailed
! formulas given in Jones and Pewsey [1] directly as published.
!
! * We could have converted the formulas to the tfp notation, but such a
! reformulation does not save significant computational time while severly
! sacrificing clarity.
!------------------------------------------------------------------------------
ELEMENTAL SUBROUTINE convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau
    REAL(8), INTENT(OUT) :: xi, eta, eps, delta

    xi    = mu
    eta   = sigma * 2.0 / SINH(ASINH(2.0) * tau)
    eps   = nu
    delta = 1.0 / tau
END SUBROUTINE


!------------------------------------------------------------------------------
! FUNCTION rational_approximation(t)
!
! Compute Abramowitz and Stegun (1965) Equation (26.2.23).
!
! Parameters
! ----------
! t : real scalar
!   The value at which to compute Equation (26.2.23).
!
! Returns
! -------
! rational_approximation : real scalar
!
! Notes
! -----
! * This is merely a helper function, for a helper function.
!
! * The rational approximation is computed using Horner's rule.
!
! * This FORTRAN code was inspired by John D. Cook's writeup found at:
! https://www.johndcook.com/blog/normal_cdf_inverse/.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION rational_approximation(t)
    REAL(8) :: rational_approximation
    REAL(8), INTENT(IN) :: t

    REAL(8), PARAMETER :: c(3) = (/ 2.515517, 0.802853, 0.010328 /)
    REAL(8), PARAMETER :: d(3) = (/ 1.432788, 0.189269, 0.001308 /)

    rational_approximation = t - ((c(3)*t + c(2))*t + c(1)) / (((d(3)*t + d(2))*t + d(1))*t + 1.0)
END FUNCTION rational_approximation


!------------------------------------------------------------------------------
! FUNCTION gaussian_cdf_inv(pr))
!
! Compute the approximate inverse cdf for a standard normal distribution.
!
! Parameters
! ----------
! pr : real scalar in the range 0 < pr < 1.
!   The value at which to compute Equation (26.2.23).
!
! Returns
! -------
! gaussian_cdf_inv : real scalar.
!   returns x where gaussian_cdf(x) = pr.
!
! Notes
! -----
! * This is merely a helper function.
!
! * This implements Abramowitz and Stegun (1965) Equation (26.2.23).
!
! * This is an approximation. The error is bounded by |error| < 4.5e-4.
!
! * If pr is out of the range 0 < pr < 1, this returns an IEEE_QUIET_NAN.
!
! * This FORTRAN code was inspired by John D. Cook's writeup found at:
! https://www.johndcook.com/blog/normal_cdf_inverse/.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION gaussian_cdf_inv(pr)
    REAL(8) :: gaussian_cdf_inv
    REAL(8), INTENT(IN) :: pr

    IF (pr <= 0.0) THEN
        gaussian_cdf_inv = -HUGE(1.0_8)
    ELSE IF (pr < 0.5) THEN
        gaussian_cdf_inv = -rational_approximation( SQRT(-2.0*LOG(pr)) )
    ELSE IF (pr < 1.0) THEN
        gaussian_cdf_inv = rational_approximation( SQRT(-2.0*LOG(1.0 - pr)) )
    ELSE
        gaussian_cdf_inv = HUGE(1.0_8)
    END IF
END FUNCTION gaussian_cdf_inv


!==============================================================================
END MODULE shash_module
!==============================================================================
