!==============================================================================
! shash.f95
!
! FORTRAN-based sinh-arcsinh normal (SHASH) distribution utility functions.
!
! Public functions
! ----------------
! pdf(x, mu, sigma, nu, tau)
!   Compute the SHASH probability density function (pdf).
!   * x, mu, sigma, nu, and tau may be scalars.
!   * x, mu, sigma, nu, and tau may be commensurate arrays.
!   * x may be an array with mu, sigma, nu, and tau scalars.
!
! cdf(x, mu, sigma, nu, tau)
!   Compute the SHASH cumulative distribution function (cdf).
!   * x, mu, sigma, nu, and tau may be scalars.
!   * x, mu, sigma, nu, and tau may be commensurate arrays.
!   * x may be an array with mu, sigma, nu, and tau scalars.
!
! quantile(pr, mu, sigma, nu, tau)
!   Compute the SHASH inverse cumulative distribution function: that is,
!   find x such that cdf(x) = pr.
!   * pr, mu, sigma, nu, and tau may be scalars.
!   * pr, mu, sigma, nu, and tau may be commensurate arrays.
!
! median(mu, sigma, nu, tau)
!   Compute the distribution median.
!   * mu, sigma, nu, and tau may be scalars.
!   * mu, sigma, nu, and tau may be commensurate arrays.
!
! mean(mu, sigma, nu, tau)
!   Compute the distribution mean.
!   * mu, sigma, nu, and tau may be scalars.
!   * mu, sigma, nu, and tau may be commensurate arrays.
!
! mode(mu, sigma, nu, tau)
!   Compute the distribution mode.
!   * mu, sigma, nu, and tau may be scalars.
!   * mu, sigma, nu, and tau may be commensurate arrays.
!
! variance(mu, sigma, nu, tau)
!   Compute the distribution variance.
!   * mu, sigma, nu, and tau may be scalars.
!   * mu, sigma, nu, and tau may be commensurate arrays.
!
! stddev(mu, sigma, nu, tau)
!   Compute the distribution standard deviation.
!   * mu, sigma, nu, and tau may be scalars.
!   * mu, sigma, nu, and tau may be commensurate arrays.
!
! skew(mu, sigma, nu, tau)
!   Compute the distribution moment coefficient of skewness.
!   * mu, sigma, nu, and tau may be scalars.
!   * mu, sigma, nu, and tau may be commensurate arrays.
!
! compute_summary_statistics(mu, sigma, nu, tau)
!   Compute a suite of summary statistics.
!   * mu, sigma, nu, and tau may be scalars.
!   * mu, sigma, nu, and tau may be commensurate arrays.
!
! Notes
! -----
! * The sinh-arcsinh normal distribution was defined in [1]. A more accessible,
!   though less comprehensive, presentation is given in [2].
!
! * This code was tested using:
!   GNU Fortran (MinGW-W64 x86_64-ucrt-posix-seh, built by Brecht Sanders) 12.2.0
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
MODULE shash_module
IMPLICIT NONE
PRIVATE

!---------------------------------------------------------------------------
! Defined types
!---------------------------------------------------------------------------
TYPE SUMMARY_STATISTICS
    REAL(8) :: interquartile_range
    REAL(8) :: mean
    REAL(8) :: median
    REAL(8) :: mode
    REAL(8) :: percentile_10
    REAL(8) :: percentile_25
    REAL(8) :: percentile_75
    REAL(8) :: percentile_90
    REAL(8) :: skew
    REAL(8) :: stddev
    REAL(8) :: variance
END TYPE

!------------------------------------------------------------------------------
! Interface block
!------------------------------------------------------------------------------
INTERFACE pdf
   MODULE PROCEDURE pdf_elemental
   MODULE PROCEDURE pdf_x_array
END INTERFACE

INTERFACE cdf
   MODULE PROCEDURE cdf_elemental
   MODULE PROCEDURE cdf_x_array
END INTERFACE

!------------------------------------------------------------------------------
! Visibility
!------------------------------------------------------------------------------
PUBLIC pdf
PUBLIC cdf
PUBLIC quantile
PUBLIC median
PUBLIC mean
PUBLIC mode
PUBLIC variance
PUBLIC stddev
PUBLIC skew

PUBLIC SUMMARY_STATISTICS
PUBLIC compute_summary_statistics

!---------------------------------------------------------------------------
! Mathematical constants
!---------------------------------------------------------------------------
REAL(8), PARAMETER :: PI                   = 3.1415926535897932384626434_8
REAL(8), PARAMETER :: PI_OVER_TWO          = 1.5707963267948966192313217_8
REAL(8), PARAMETER :: ONE_OVER_SQRT_TWO    = 0.7071067811865475244008444_8
REAL(8), PARAMETER :: ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399461_8


CONTAINS

!------------------------------------------------------------------------------
! FUNCTION pdf_elemental(x, mu, sigma, nu, tau)
!
! Compute the SHASH probability density function (pdf).
!
! Arguments
! ---------
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
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION pdf_elemental(x, mu, sigma, nu, tau) RESULT(pdf)
    REAL(8) :: pdf
    REAL(8), INTENT(IN) :: x, mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: y, z

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) Equation (2) on page 762.
    y = (x - xi) / eta
    z = SINH(delta * ASINH(y) - eps)
    pdf = ONE_OVER_SQRT_TWO_PI * (delta / eta) * SQRT((1.0 + z*z) / (1.0 + y*y)) * EXP(-z*z / 2.0)
END FUNCTION pdf_elemental


!------------------------------------------------------------------------------
! FUNCTION pdf_x_array(x, mu, sigma, nu, tau)
!
! Compute the SHASH probability density function (pdf) for an array x.
!
! Arguments
! ---------
! x : real array
!   The values at which to compute the SHASH probability density function.
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
! pdf : real array
!   The computed shash(mu, sigma, nu, tau) probability density function
!   evaluated at x.  pdf is the same size as x.
!------------------------------------------------------------------------------
PURE FUNCTION pdf_x_array(x, mu, sigma, nu, tau) RESULT(pdf)
    REAL(8), INTENT(IN) :: x(:)
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: pdf

    REAL(8) :: xi, eta, eps, delta
    REAL(8), DIMENSION( SIZE(x) ) :: y, z

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) Equation (2) on page 762.
    y = (x - xi) / eta
    z = SINH(delta * ASINH(y) - eps)
    pdf = ONE_OVER_SQRT_TWO_PI * (delta / eta) * SQRT((1.0 + z*z) / (1.0 + y*y)) * EXP(-z*z / 2.0)
END FUNCTION pdf_x_array


!------------------------------------------------------------------------------
! FUNCTION cdf_elemental(x, mu, sigma, nu, tau)
!
! Compute the SHASH cumulative distribution function (cdf).
!
! Arguments
! ---------
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
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION cdf_elemental(x, mu, sigma, nu, tau) RESULT(CDF)
    REAL(8) :: cdf
    REAL(8), INTENT(IN) :: x, mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: y, z

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762.
    y = (x - xi) / eta
    z = SINH(delta * ASINH(y) - eps)
    cdf = 0.5 * (1.0 + ERF(ONE_OVER_SQRT_TWO * z))
END FUNCTION cdf_elemental


!------------------------------------------------------------------------------
! FUNCTION cdf_x_array(x, mu, sigma, nu, tau)
!
! Compute the SHASH cumulative distribution function (cdf).
!
! Arguments
! ---------
! x : real array
!   The values at which to compute the SHASH cumulative distribution function.
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
! cdf : real array
!   The computed shash(mu, sigma, nu, tau) cumulative distribution function
!   evaluated at x. cdf is the same size as x.
!------------------------------------------------------------------------------
PURE FUNCTION cdf_x_array(x, mu, sigma, nu, tau) RESULT(cdf)
    REAL(8), INTENT(IN) :: x(:)
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau
    REAL(8), DIMENSION( SIZE(x) ) :: cdf

    REAL(8) :: xi, eta, eps, delta
    REAL(8), DIMENSION( SIZE(x) ) :: y, z

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762.
    y = (x - xi) / eta
    z = SINH(delta * ASINH(y) - eps)
    cdf = 0.5 * (1.0 + ERF(ONE_OVER_SQRT_TWO * z))
END FUNCTION cdf_x_array


!------------------------------------------------------------------------------
! FUNCTION quantile(pr, mu, sigma, nu, tau)
!
! Compute the SHASH inverse cumulative distribution function: that is,
! find x such that cdf(x) = pr.
!
! Arguments
! ---------
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
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
! * This function uses Newton's method, with a pretty good starting guess.
!   The result is accurate (|error in pr| < 1e-6), but it is a bit slow.
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
! FUNCTION mean(mu, sigma, nu, tau)
!
! Compute the SHASH distribution mean.
!
! Arguments
! ---------
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
! mean : real scalar
!   The computed mean of shash(mu, sigma, nu, tau).
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION mean(mu, sigma, nu, tau)
    REAL(8) :: mean
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: evX

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) middle of page 764.
    evX = SINH(eps / delta) * jones_pewsey_p(1.0 / delta)
    mean = xi + eta * evX
END FUNCTION mean


!------------------------------------------------------------------------------
! FUNCTION median(mu, sigma, nu, tau)
!
! Compute the SHASH distribution median.
!
! Arguments
! ---------
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
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
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


!------------------------------------------------------------------------------
! FUNCTION mode(mu, sigma, nu, tau)
!
! Compute the SHASH distribution mode.
!
! Arguments
! ---------
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
! mode : real scalar
!   The computed mode of shash(mu, sigma, nu, tau).
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
! * We cannot solve for the mode analytically, so we compute it
!   numerically using Newton's method to find the root of the derivative
!   of the density function.
!
! * We initialize the Newton's method iterations using the medians.
!
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION mode(mu, sigma, nu, tau)
    REAL(8) :: mode
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: y, S, C, T, g, dg, deltay
    INTEGER :: iteration

    INTEGER, PARAMETER :: MAX_ITERATION = 100
    REAL(8), PARAMETER :: TOLERANCE = 1.0e-4

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Use the median as the initial guess.
    y = SINH(eps / delta)

    ! Apply Newton's method.
    DO iteration = 1, MAX_ITERATION
        S = SINH(delta * ASINH(y) - eps)
        C = COSH(delta * ASINH(y) - eps)
        T = SQRT(1.0 + y**2)

        ! g is the derivative of the pdf. dg is the derivative of g.
        g = y + y * S**2 + delta * C * S**3 * T
        dg = 1.0 + S**2 + delta**2 * S**2 * (3.0*C**2 + S**2) + 2.0*delta*C*S*(1.0 + S**2) * y/T

        deltay = g/dg
        y = y - deltay
        IF (ABS(deltay) < TOLERANCE) EXIT
    END DO

    mode = xi + eta * y
END FUNCTION mode


!------------------------------------------------------------------------------
! FUNCTION variance(mu, sigma, nu, tau)
!
! Compute the SHASH distribution variance.
!
! Arguments
! ---------
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
! variance : real scalar
!   The computed variance of shash(mu, sigma, nu, tau).
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
! * This code uses two basic formulas:
!
!   var(X) = E(X^2) - (E(X))^2
!   var(a*X + b) = a^2 * var(X)
!
! * The E(X) and E(X^2) are computed using the moment equations given on
!   page 764 of [1].
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION variance(mu, sigma, nu, tau)
    REAL(8) :: variance
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: evX, evX2

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) middle of page 764.
    evX = SINH(eps / delta) * jones_pewsey_p(1.0 / delta)
    evX2 = (COSH(2.0 * eps / delta) * jones_pewsey_p(2.0 / delta) - 1.0) / 2.0
    variance = eta**2 * (evX2 - evX**2)
END FUNCTION variance


!------------------------------------------------------------------------------
! FUNCTION stddev(mu, sigma, nu, tau)
!
! Compute the SHASH distribution standard deviation.
!
! Arguments
! ---------
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
! stddev : real scalar
!   The computed standard deviation of shash(mu, sigma, nu, tau).
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION stddev(mu, sigma, nu, tau)
    REAL(8) :: stddev
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau

    stddev = SQRT(variance(mu, sigma, nu, tau))
END FUNCTION stddev


!------------------------------------------------------------------------------
! FUNCTION skew(mu, sigma, nu, tau)
!
! Compute the SHASH distribution moment coefficient of skewness.
!
! Arguments
! ---------
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
! skew : real scalar
!   The computed skewness of shash(mu, sigma, nu, tau).
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
! * This function computes the "moment coefficient of skewness", also known as
!   the "third standardized moment", also known as "Fisher's moment coefficient
!   of skewness", also known as "Pearson's moment coefficient of skewness".
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION skew(mu, sigma, nu, tau)
    REAL(8) :: skew
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau

    REAL(8) :: xi, eta, eps, delta
    REAL(8) :: evX, evX2, evX3
    REAL(8) :: evY, evY3, stdY

    CALL convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) middle of page 764.
    evX = SINH(eps / delta) * jones_pewsey_p(1.0 / delta)
    evX2 = (COSH(2.0 * eps / delta) * jones_pewsey_p(2.0 / delta) - 1.0) / 2.0
    evX3 = (SINH(3.0 * eps / delta) * jones_pewsey_p(3.0 / delta) - 3.0 * SINH(eps / delta) * jones_pewsey_p(1.0 / delta)) / 4.0

    evY = xi + eta * evX
    evY3 = xi**3 + 3.0 * xi**2 * eta * evX + 3.0 * xi * eta**2 * evX2 + eta**3 * evX3
    stdY = eta * SQRT(evX2 - evX**2)

    skew = (evY3 - 3.0 * evY * stdY**2 - evY**3) / stdY**3
END FUNCTION skew


!------------------------------------------------------------------------------
! ELEMENTAL FUNCTION compute_summary_statistics(mu, sigma, nu, tau, summary)
!
! Compute a set of summary statistics for the SHASH distribution
!
! Arguments
! ---------
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
! summary : SUMMARY_STATISTICS, on output
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION compute_summary_statistics(mu, sigma, nu, tau) RESULT(summary)
    REAL(8), INTENT(IN) :: mu, sigma, nu, tau
    TYPE (SUMMARY_STATISTICS) :: summary

    summary%mean   = mean(mu, sigma, nu, tau)
    summary%median = median(mu, sigma, nu, tau)
    summary%mode   = mode(mu, sigma, nu, tau)

    summary%percentile_10 = quantile(0.10_8, mu, sigma, nu, tau)
    summary%percentile_25 = quantile(0.25_8, mu, sigma, nu, tau)
    summary%percentile_75 = quantile(0.75_8, mu, sigma, nu, tau)
    summary%percentile_90 = quantile(0.90_8, mu, sigma, nu, tau)

    summary%interquartile_range = summary%percentile_75 - summary%percentile_25
    summary%skew = skew(mu, sigma, nu, tau)

    summary%stddev = stddev(mu, sigma, nu, tau)
    summary%variance = variance(mu, sigma, nu, tau)
END FUNCTION compute_summary_statistics


!==============================================================================
! PRIVATE HELPER FUNCTIONS
!==============================================================================

!------------------------------------------------------------------------------
! SUBROUTINE convert_tf_to_jp(mu, sigma, nu, tau, xi, eta, eps, delta)
!
! Convert from TensorFlow formulation and notation to the Jones and Pewsey
! formulation and notation.
!
! Arguments
! ---------
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
!   formulas given in Jones and Pewsey [1] directly as published.
!
! * We could have converted the formulas to the tfp notation, but such a
!   reformulation does not save significant computational time while severly
!   sacrificing clarity.
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
! FUNCTION gaussian_cdf_inv(pr))
!
! Compute the approximate inverse cdf for a standard normal distribution.
!
! Arguments
! ---------
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
! * If pr is out of the range 0 < pr < 1, this routine fails.
!
! * This FORTRAN code was inspired by John D. Cook's writeup found at:
!   https://www.johndcook.com/blog/normal_cdf_inverse/.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION gaussian_cdf_inv(pr) RESULT(cdf_inv)
    REAL(8) :: cdf_inv
    REAL(8), INTENT(IN) :: pr

    IF (pr <= 0.0) THEN
        cdf_inv = -HUGE(1.0_8)
    ELSE IF (pr < 0.5) THEN
        cdf_inv = -rational_approximation( SQRT(-2.0*LOG(pr)) )
    ELSE IF (pr < 1.0) THEN
        cdf_inv = rational_approximation( SQRT(-2.0*LOG(1.0 - pr)) )
    ELSE
        cdf_inv = HUGE(1.0_8)
    END IF
END FUNCTION gaussian_cdf_inv


!------------------------------------------------------------------------------
! FUNCTION rational_approximation(t)
!
! Compute an approximate inverse cdf for a standard normal distribution using
! the approximation given by Abramowitz and Stegun (1965) Equation (26.2.23).
!
! Arguments
! ---------
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
!   https://www.johndcook.com/blog/normal_cdf_inverse/.
!
! References
! ----------
! Abramowitz, M. & Stegun, I. A., Handbook of Mathematical Functions:
! with Formulas, Graphs, and Mathematical Tables, Dover Publications, 1965.
! ISBN: 978-0486612720.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION rational_approximation(t)
    REAL(8) :: rational_approximation
    REAL(8), INTENT(IN) :: t

    REAL(8), PARAMETER :: c(3) = (/ 2.515517, 0.802853, 0.010328 /)
    REAL(8), PARAMETER :: d(3) = (/ 1.432788, 0.189269, 0.001308 /)

    rational_approximation = t - ((c(3)*t + c(2))*t + c(1)) / (((d(3)*t + d(2))*t + d(1))*t + 1.0)
END FUNCTION rational_approximation


!------------------------------------------------------------------------------
! FUNCTION jones_pewsey_p(q)
!
! Compute the Jones and Pewsey factor Pq as defined on page 764 of Jones and
! Pewsey (2009).
!
! Arguments
! ---------
! q : real scalar, q > 0
!
! Returns
! -------
! p : real scalar
!   Jones and Pewsey Pq factor.
!
! Notes
! -----
! * The Pq factor is used in the computation of the moments for a
!   sinh-arcsinh normal (SHASH) random variable.
!
! * The Jones and Pewsey Pq factor is defined on page 764 of [1].
!
!   $P_q = \frac{e^{1/4}}{(8 \pi)^{1/2}} \cdot
!          \left[ K_{(q+1)/2}(1/4) + K_{(q-1)/2}(1/4) \right]$
!
!   where $K_{v}(1/4)$ is the modified Bessel function of the second kind.
!
! * The strange constant 0.25612..., is $\frac{e^{1/4}}{(8\pi)^{1/2}}$
!   computed using a high-precision calculator.
!
! References
! ----------
! [1] M. C. Jones and A. Pewsey. "Sinh-arcsinh distributions", Biometrika 96.4
! October 2009, pp. 761–780. DOI: 10.1093/biomet/asp053.
!
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION jones_pewsey_p(q) RESULT(f)
    REAL(8) :: f
    REAL(8), INTENT(IN) :: q

    f = 0.2561260139134036986353746_8 * (Kv((q+1.0)/2.0) + Kv((q-1.0)/2.0))
END FUNCTION jones_pewsey_p


!------------------------------------------------------------------------------
! FUNCTION Kv(v)
!
! Compute the modified Bessel function of the second kind at (1/4):
! $K_{v}(1/4)$.
!
! Arguments
! ---------
! v :: real scalar, v > -1
!
! Returns
! -------
! besselk : real scalar
!   The value of $K_{v}(1/4)$.
!
! Notes
! -----
! * For all but near-integer values of v, the computation uses Equation
!   (6.1.4) of page 203 of [1].
!
!       $K_v = \frac{\pi}{2} \frac{I_{-v} - I_{v}}{\sin(v \pi)}$
!
! * The standard relationship between Kv(z) and Iv(z) has computational
!   singularities at v = 0, 1, 2, ..., due to the sin(v*PI) term in the
!   denominator. To eliminate these computational singularities we use
!   a quadratic approximation around v = 0, a linear approximation around
!   v = 1, and the standard forward recursion around v = 2, 3, ...
!
!       $K_{v-1} + 8v K_{v} = K_{v+1}$
!
!   See Equation (61.23) on page 206 of [1], and substitute in z = 1/4.
!
! * The quadratic approximation about v=0, and the linear approximation
!   around v=1, are the truncated Taylor series.
!
! References
! ----------
! [1] S. Zhang and J. Jin. Computation of Special Functions. John Wiley and
! Sons, Inc., 1996. ISBN 978-0471119630.
!
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION Kv(v) RESULT(besselk)
    REAL(8) :: besselk
    REAL(8), INTENT(IN) :: v

    REAL(8), PARAMETER :: TOLERANCE = 1e-4
    REAL(8) :: nu, besselk_minus1, besselk_plus1
    INTEGER :: j

    IF (ABS(v - NINT(v)) > TOLERANCE) THEN
        besselk = PI_OVER_TWO * (Iv(-v) - Iv(v)) / SIN(v*PI)
    ELSEIF (ABS(v) < TOLERANCE) THEN
        besselk = 1.541506751248303 + 1.491844425771660 * v**2
    ELSEIF (ABS(v - 1.0) < TOLERANCE) THEN
        besselk = 3.747025974440712 + 6.166027005560792 * (v-1.0)
    ELSE
        nu = 1.0 + v - NINT(v)

        besselk_minus1 = 1.541506751248303 + 1.491844425771660 * (v - NINT(v))**2
        besselk        = 3.747025974440712 + 6.166027005560792 * (v - NINT(v))

        DO j = 1, (NINT(v) - 1)
            besselk_plus1 = besselk_minus1 + 8.0*nu*besselk

            nu = nu + 1.0
            besselk_minus1 = besselk
            besselk = besselk_plus1
        END DO
    END IF
END FUNCTION Kv


!------------------------------------------------------------------------------
! FUNCTION Iv(v)
!
! Compute the modified Bessel function of the first kind at (1/4):
! $I_{v}(1/4)$.
!
! Arguments
! ---------
! v :: real scalar, v > -1
!
! Returns
! -------
! f : real scalar
!   The value of $I_{v}(1/4)$.
!
! Notes
! -----
! * The computation is carried out using
!
!       $I_v = \sum_{m=0}^{\infty} \frac{1}{8^{2m+v} m! \Gamma(m+v+1)}
!
!   See Equation (6.1.2) on page 202 of [1], and substitute in z = 1/4.
!
! * Using only six terms in the infinite series is sufficient to achieve full
!   double precision for all necessary cases; i.e. for $-1/2 < v < \infty$.
!
! References
! ----------
! [1] S. Zhang and J. Jin. Computation of Special Functions. John Wiley and
! Sons, Inc., 1996. ISBN 978-0471119630.
!
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION Iv(v) RESULT(f)
    REAL(8) :: f
    REAL(8), INTENT(IN) :: v

    REAL(8) :: term
    INTEGER :: m

    term = 1.0 / (8.0**v * GAMMA(v + 1.0))
    f = term

    DO m = 1, 6
        term = term / (64.0 * m * (v+m))
        f = f + term
    END DO
END FUNCTION Iv


!==============================================================================
END MODULE shash_module
!==============================================================================
