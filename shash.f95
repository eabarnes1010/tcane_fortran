!==============================================================================
! shash.f95
!
! FORTRAN-based sinh-arcsinh normal (SHASH) distribution utility functions.
!
! Public functions
! ----------------
! pdf(x, loc, sigma, skewness, tailweight)
!   Compute the SHASH probability density function (pdf).
!   * x, loc, sigma, skewness, and tailweight may be scalars.
!   * x, loc, sigma, skewness, and tailweight may be commensurate arrays.
!   * x may be an array with loc, sigma, skewness, and tailweight scalars.
!
! cdf(x, loc, sigma, skewness, tailweight)
!   Compute the SHASH cumulative distribution function (cdf).
!   * x, loc, sigma, skewness, and tailweight may be scalars.
!   * x, loc, sigma, skewness, and tailweight may be commensurate arrays.
!   * x may be an array with loc, sigma, skewness, and tailweight scalars.
!
! quantile(pr, loc, sigma, skewness, tailweight)
!   Compute the SHASH inverse cumulative distribution function: that is,
!   find x such that cdf(x) = pr.
!   * pr, loc, sigma, skewness, and tailweight may be scalars.
!   * pr, loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! median(loc, sigma, skewness, tailweight)
!   Compute the distribution median.
!   * loc, sigma, skewness, and tailweight may be scalars.
!   * loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! mean(loc, sigma, skewness, tailweight)
!   Compute the distribution mean.
!   * loc, sigma, skewness, and tailweight may be scalars.
!   * loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! mode(loc, sigma, skewness, tailweight)
!   Compute the distribution mode.
!   * loc, sigma, skewness, and tailweight may be scalars.
!   * loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! variance(loc, sigma, skewness, tailweight)
!   Compute the distribution variance.
!   * loc, sigma, skewness, and tailweight may be scalars.
!   * loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! stddev(loc, sigma, skewness, tailweight)
!   Compute the distribution standard deviation.
!   * loc, sigma, skewness, and tailweight may be scalars.
!   * loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! skew(loc, sigma, skewness, tailweight)
!   Compute the distribution moment coefficient of skewness.
!   * loc, sigma, skewness, and tailweight may be scalars.
!   * loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! compute_summary_statistics(loc, sigma, skewness, tailweight)
!   Compute a suite of summary statistics.
!   * loc, sigma, skewness, and tailweight may be scalars.
!   * loc, sigma, skewness, and tailweight may be commensurate arrays.
!
! Notes
! -----
! * The sinh-arcsinh normal distribution was defined in [1]. A more
!   accessible, though less comprehensive, presentation is given in [2].
!
! * The argument "sigma" is the same as the Tensorflow "scale". We could not
!   use "scale" because "scale" is the name of an intrinsic FORTRAN function.
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
module shash_module
implicit none
private

!------------------------------------------------------------------------------
! Defined types
!------------------------------------------------------------------------------
type summary_statistics
    real(8) :: interquartile_range
    real(8) :: mean
    real(8) :: median
    real(8) :: mode
    real(8) :: percentile_10
    real(8) :: percentile_25
    real(8) :: percentile_75
    real(8) :: percentile_90
    real(8) :: skew
    real(8) :: stddev
    real(8) :: variance
end type

!------------------------------------------------------------------------------
! Interface block
!------------------------------------------------------------------------------
interface pdf
   module procedure pdf_elemental
   module procedure pdf_x_array
end interface

interface cdf
   module procedure cdf_elemental
   module procedure cdf_x_array
end interface

!------------------------------------------------------------------------------
! Visibility
!------------------------------------------------------------------------------
public pdf
public cdf
public quantile
public median
public mean
public mode
public variance
public stddev
public skew

public summary_statistics
public compute_summary_statistics

!------------------------------------------------------------------------------
! Mathematical constants
!------------------------------------------------------------------------------
real(8), parameter :: PI                   = 3.1415926535897932384626434_8
real(8), parameter :: PI_OVER_TWO          = 1.5707963267948966192313217_8
real(8), parameter :: ONE_OVER_SQRT_TWO    = 0.7071067811865475244008444_8
real(8), parameter :: ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399461_8


contains

!------------------------------------------------------------------------------
! function pdf_elemental(x, loc, sigma, skewness, tailweight)
!
! Compute the SHASH probability density function (pdf).
!
! Arguments
! ---------
! x : real scalar
!   The value at which to compute the SHASH probability density function.
!
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! pdf : real scalar
!   The computed probability density function evaluated at x.
!
!------------------------------------------------------------------------------
elemental function pdf_elemental(x, loc, sigma, skewness, tailweight) result(pdf)
    real(8) :: pdf
    real(8), intent(in) :: x, loc, sigma, skewness, tailweight

    real(8) :: xi, eta, eps, delta
    real(8) :: y, z

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) Equation (2) on page 762.
    y = (x - xi) / eta
    z = sinh(delta * asinh(y) - eps)
    pdf = ONE_OVER_SQRT_TWO_PI * (delta / eta) * sqrt((1.0 + z**2) / (1.0 + y**2)) * exp(-z**2 / 2.0)
end function pdf_elemental


!------------------------------------------------------------------------------
! function pdf_x_array(x, loc, sigma, skewness, tailweight)
!
! Compute the SHASH probability density function (pdf) for an array x.
!
! Arguments
! ---------
! x : real array
!   The values at which to compute the SHASH probability density function.
!
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! pdf : real array
!   The computed probability density function evaluated at x.  The returned
!   pdf is the same size as x.
!
!------------------------------------------------------------------------------
pure function pdf_x_array(x, loc, sigma, skewness, tailweight) result(pdf)
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(x) ) :: pdf

    real(8) :: xi, eta, eps, delta
    real(8), dimension( size(x) ) :: y, z

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) Equation (2) on page 762.
    y = (x - xi) / eta
    z = sinh(delta * asinh(y) - eps)
    pdf = ONE_OVER_SQRT_TWO_PI * (delta / eta) * sqrt((1.0 + z*z) / (1.0 + y*y)) * exp(-z*z / 2.0)
end function pdf_x_array


!------------------------------------------------------------------------------
! function cdf_elemental(x, loc, sigma, skewness, tailweight)
!
! Compute the SHASH cumulative distribution function (cdf).
!
! Arguments
! ---------
! x : real scalar
!   The value at which to compute the SHASH cumulative distribution function.
!
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! cdf : real scalar
!   The computed cumulative distribution function evaluated at x.
!
!------------------------------------------------------------------------------
elemental function cdf_elemental(x, loc, sigma, skewness, tailweight) result(cdf)
    real(8), intent(in) :: x, loc, sigma, skewness, tailweight
    real(8) :: cdf

    real(8) :: xi, eta, eps, delta
    real(8) :: y, z

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762.
    y = (x - xi) / eta
    z = sinh(delta * asinh(y) - eps)
    cdf = 0.5 * (1.0 + erf(ONE_OVER_SQRT_TWO * z))
end function cdf_elemental


!------------------------------------------------------------------------------
! function cdf_x_array(x, loc, sigma, skewness, tailweight)
!
! Compute the SHASH cumulative distribution function (cdf).
!
! Arguments
! ---------
! x : real array
!   The values at which to compute the SHASH cumulative distribution function.
!
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! cdf : real array
!   The computed cumulative distribution function evaluated at x. The returned
!   cdf is the same size as x.
!
!------------------------------------------------------------------------------
pure function cdf_x_array(x, loc, sigma, skewness, tailweight) result(cdf)
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(x) ) :: cdf

    real(8) :: xi, eta, eps, delta
    real(8), dimension( size(x) ) :: y, z

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762.
    y = (x - xi) / eta
    z = sinh(delta * asinh(y) - eps)
    cdf = 0.5 * (1.0 + erf(ONE_OVER_SQRT_TWO * z))
end function cdf_x_array


!------------------------------------------------------------------------------
! function quantile(pr, loc, sigma, skewness, tailweight)
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
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! quantile : real scalar
!   The computed inverse cumulative distribution function evaluated
!   at probability pr.
!
! Notes
! -----
! * This function uses Newton's method, with a pretty good starting guess.
!   The result is accurate (|error in pr| < 1e-6), but it is a bit slow.
!
! * If pr is out of the range 0 < pr < 1, this function fails.
!------------------------------------------------------------------------------
elemental function quantile(pr, loc, sigma, skewness, tailweight) result(x)
    real(8), intent(in) :: pr, loc, sigma, skewness, tailweight
    real(8) :: x

    real(8) :: xi, eta, eps, delta, z

    integer :: iteration
    integer, parameter :: MAX_ITERATION = 5

    real(8) :: difference
    real(8), parameter :: MAX_ABS_DIFFERENCE = 1.0e-6

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762, and an approximate
    ! Gaussian cdf inverse, to get an initial approximation.
    z = gaussian_cdf_inv(pr)
    x = xi + eta * sinh((asinh(z) + eps) / delta)

    ! Improve the approximation using Newton's method.
    do iteration = 1, MAX_ITERATION
        difference = cdf(x, loc, sigma, skewness, tailweight) - pr
        x = x - difference/pdf(x, loc, sigma, skewness, tailweight)
        if (abs(difference) <= MAX_ABS_DIFFERENCE) exit
    end do
end function quantile


!------------------------------------------------------------------------------
! function mean(loc, sigma, skewness, tailweight)
!
! Compute the SHASH distribution mean.
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! mean : real scalar
!   The computed mean.
!
!------------------------------------------------------------------------------
elemental function mean(loc, sigma, skewness, tailweight)
    real(8) :: mean
    real(8), intent(in) :: loc, sigma, skewness, tailweight

    real(8) :: xi, eta, eps, delta
    real(8) :: evX

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) middle of page 764.
    evX = sinh(eps / delta) * jones_pewsey_p(1.0 / delta)
    mean = xi + eta * evX
end function mean


!------------------------------------------------------------------------------
! function median(loc, sigma, skewness, tailweight)
!
! Compute the SHASH distribution median.
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! median : real scalar
!   The computed median.
!
! Notes
! -----
! * Unlike the quantile function, this is an exact computation. This
!   computation is also much faster than the quantile function.  As such, this
!   function should be used in lieu of quantile(0.5, ...).
!
!------------------------------------------------------------------------------
elemental function median(loc, sigma, skewness, tailweight)
    real(8) :: median
    real(8), intent(in) :: loc, sigma, skewness, tailweight

    real(8) :: xi, eta, eps, delta

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) bottom of page 762.
    median = xi + eta * sinh(eps / delta)
end function median


!------------------------------------------------------------------------------
! function mode(loc, sigma, skewness, tailweight)
!
! Compute the SHASH distribution mode.
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! mode : real scalar
!   The computed mode.
!
! Notes
! -----
! * We cannot solve for the mode analytically, so we compute it
!   numerically using Newton's method to find the root of the derivative
!   of the density function.
!
! * We initialize the Newton's method iterations using the median.
!
!------------------------------------------------------------------------------
elemental function mode(loc, sigma, skewness, tailweight)
    real(8) :: mode
    real(8), intent(in) :: loc, sigma, skewness, tailweight

    real(8) :: xi, eta, eps, delta
    real(8) :: y, S, C, T, g, dg, deltay
    integer :: iteration

    integer, parameter :: MAX_ITERATION = 100
    real(8), parameter :: TOLERANCE = 1.0e-4

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Use the median as the initial guess.
    y = sinh(eps / delta)

    ! Apply Newton's method.
    do iteration = 1, MAX_ITERATION
        S = sinh(delta * asinh(y) - eps)
        C = cosh(delta * asinh(y) - eps)
        T = sqrt(1.0 + y**2)

        ! g is the derivative of the pdf. dg is the derivative of g.
        g = y + y * S**2 + delta * C * S**3 * T
        dg = 1.0 + S**2 + delta**2 * S**2 * (3.0*C**2 + S**2) + 2.0*delta*C*S*(1.0 + S**2) * y/T

        deltay = g/dg
        y = y - deltay
        if (abs(deltay) < TOLERANCE) exit
    end do

    mode = xi + eta * y
end function mode


!------------------------------------------------------------------------------
! function variance(loc, sigma, skewness, tailweight)
!
! Compute the SHASH distribution variance.
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! variance : real scalar
!   The computed variance.
!
! Notes
! -----
! * This code uses two basic formulas:
!
!   var(X) = E(X^2) - (E(X))^2
!   var(a*X + b) = a^2 * var(X)
!
! * The E(X) and E(X^2) are computed using the moment equations given on
!   page 764 of [1].
!
!------------------------------------------------------------------------------
elemental function variance(loc, sigma, skewness, tailweight)
    real(8) :: variance
    real(8), intent(in) :: loc, sigma, skewness, tailweight

    real(8) :: xi, eta, eps, delta
    real(8) :: evX, evX2

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) middle of page 764.
    evX = sinh(eps / delta) * jones_pewsey_p(1.0 / delta)
    evX2 = (cosh(2.0 * eps / delta) * jones_pewsey_p(2.0 / delta) - 1.0) / 2.0
    variance = eta**2 * (evX2 - evX**2)
end function variance


!------------------------------------------------------------------------------
! function stddev(loc, sigma, skewness, tailweight)
!
! Compute the SHASH distribution standard deviation.
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! stddev : real scalar
!   The computed standard deviation.
!
!------------------------------------------------------------------------------
elemental function stddev(loc, sigma, skewness, tailweight)
    real(8) :: stddev
    real(8), intent(in) :: loc, sigma, skewness, tailweight

    stddev = sqrt(variance(loc, sigma, skewness, tailweight))
end function stddev


!------------------------------------------------------------------------------
! function skew(loc, sigma, skewness, tailweight)
!
! Compute the SHASH distribution moment coefficient of skewness.
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! skew : real scalar
!   The computed moment coefficient of skewness.
!
! Notes
! -----
! * The moment coefficient of skewness is the third standardized moment. It
!   is also known as "Pearson's moment coefficient of skewness".
!
!------------------------------------------------------------------------------
elemental function skew(loc, sigma, skewness, tailweight)
    real(8) :: skew
    real(8), intent(in) :: loc, sigma, skewness, tailweight

    real(8) :: xi, eta, eps, delta
    real(8) :: evX, evX2, evX3
    real(8) :: evY, evY3, stdY

    call convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)

    ! Apply Jones and Pewsey (2009) middle of page 764.
    evX = sinh(eps / delta) * jones_pewsey_p(1.0 / delta)
    evX2 = (cosh(2.0 * eps / delta) * jones_pewsey_p(2.0 / delta) - 1.0) / 2.0
    evX3 = (sinh(3.0 * eps / delta) * jones_pewsey_p(3.0 / delta) - 3.0 * sinh(eps / delta) * jones_pewsey_p(1.0 / delta)) / 4.0

    evY = xi + eta * evX
    evY3 = xi**3 + 3.0 * xi**2 * eta * evX + 3.0 * xi * eta**2 * evX2 + eta**3 * evX3
    stdY = eta * sqrt(evX2 - evX**2)

    skew = (evY3 - 3.0 * evY * stdY**2 - evY**3) / stdY**3
end function skew


!------------------------------------------------------------------------------
! ELEMENTAL function compute_summary_statistics(loc, sigma, skewness, tailweight, summary)
!
! Compute a set of summary statistics for the SHASH distribution
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Returns
! -------
! summary : summary_statistics, on output
!
!------------------------------------------------------------------------------
elemental function compute_summary_statistics(loc, sigma, skewness, tailweight) result(summary)
    real(8), intent(in) :: loc, sigma, skewness, tailweight
    type (summary_statistics) :: summary

    summary%mean   = mean(loc, sigma, skewness, tailweight)
    summary%median = median(loc, sigma, skewness, tailweight)
    summary%mode   = mode(loc, sigma, skewness, tailweight)

    summary%percentile_10 = quantile(0.10_8, loc, sigma, skewness, tailweight)
    summary%percentile_25 = quantile(0.25_8, loc, sigma, skewness, tailweight)
    summary%percentile_75 = quantile(0.75_8, loc, sigma, skewness, tailweight)
    summary%percentile_90 = quantile(0.90_8, loc, sigma, skewness, tailweight)

    summary%interquartile_range = summary%percentile_75 - summary%percentile_25
    summary%skew = skew(loc, sigma, skewness, tailweight)

    summary%stddev = stddev(loc, sigma, skewness, tailweight)
    summary%variance = variance(loc, sigma, skewness, tailweight)
end function compute_summary_statistics


!==============================================================================
! PRIVATE HELPER FUNCTIONS
!==============================================================================

!------------------------------------------------------------------------------
! SUBROUTINE convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)
!
! Convert from TensorFlow formulation and notation to the Jones and Pewsey
! formulation and notation.
!
! Arguments
! ---------
! loc : real scalar
!   The location parameter.
!
! sigma : real scalar, sigma > 0
!   The scale parameter.
!
! skewness : real scalar
!   The skewness parameter.
!
! tailweight : real scalar, tailweight > 0
!   The tail-weight parameter.
!
! Notes
! -----
! * This conversion of the parameters is necessary to use the detailed
!   formulas given in Jones and Pewsey [1] directly as published.
!
! * We could have converted the formulas to the tfp notation, but such a
!   reformulation does not save significant computational time while severely
!   sacrificing clarity.
!
!------------------------------------------------------------------------------
elemental subroutine convert_tf_to_jp(loc, sigma, skewness, tailweight, xi, eta, eps, delta)
    real(8), intent(in) :: loc, sigma, skewness, tailweight
    real(8), intent(out) :: xi, eta, eps, delta

    xi    = loc
    eta   = sigma * 2.0 / sinh(asinh(2.0) * tailweight)
    eps   = skewness
    delta = 1.0 / tailweight
end subroutine


!------------------------------------------------------------------------------
! function gaussian_cdf_inv(pr))
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
! * This implements Abramowitz and Stegun (1965) Equation (26.2.23).
!
! * This is an approximation. The error is bounded by |error| < 4.5e-4.
!
! * If pr is out of the range 0 < pr < 1, this routine fails.
!
! * This FORTRAN code was inspired by John D. Cook's writeup found at:
!   https://www.johndcook.com/blog/normal_cdf_inverse/.
!
!------------------------------------------------------------------------------
elemental function gaussian_cdf_inv(pr) result(cdf_inv)
    real(8) :: cdf_inv
    real(8), intent(in) :: pr

    if (pr <= 0.0) then
        cdf_inv = -huge(1.0_8)
    else if (pr < 0.5) then
        cdf_inv = -rational_approximation( sqrt(-2.0*log(pr)) )
    else if (pr < 1.0) then
        cdf_inv = rational_approximation( sqrt(-2.0*log(1.0 - pr)) )
    else
        cdf_inv = huge(1.0_8)
    end if
end function gaussian_cdf_inv


!------------------------------------------------------------------------------
! function rational_approximation(t)
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
!
!------------------------------------------------------------------------------
elemental function rational_approximation(t)
    real(8) :: rational_approximation
    real(8), intent(in) :: t

    real(8), parameter :: c(3) = [2.515517, 0.802853, 0.010328]
    real(8), parameter :: d(3) = [1.432788, 0.189269, 0.001308]

    rational_approximation = t - ((c(3)*t + c(2))*t + c(1)) / (((d(3)*t + d(2))*t + d(1))*t + 1.0)
end function rational_approximation


!------------------------------------------------------------------------------
! function jones_pewsey_p(q)
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
elemental function jones_pewsey_p(q) result(f)
    real(8) :: f
    real(8), intent(in) :: q

    f = 0.2561260139134036986353746_8 * (Kv((q+1.0)/2.0) + Kv((q-1.0)/2.0))
end function jones_pewsey_p


!------------------------------------------------------------------------------
! function Kv(v)
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
elemental function Kv(v) result(besselk)
    real(8) :: besselk
    real(8), intent(in) :: v

    real(8), parameter :: TOLERANCE = 1e-4
    real(8) :: skewness, besselk_minus1, besselk_plus1
    integer :: j

    if (abs(v - nint(v)) > TOLERANCE) then
        besselk = PI_OVER_TWO * (Iv(-v) - Iv(v)) / sin(v*PI)
    elseif (abs(v) < TOLERANCE) then
        besselk = 1.541506751248303 + 1.491844425771660 * v**2
    elseif (abs(v - 1.0) < TOLERANCE) then
        besselk = 3.747025974440712 + 6.166027005560792 * (v-1.0)
    else
        skewness = 1.0 + v - nint(v)

        besselk_minus1 = 1.541506751248303 + 1.491844425771660 * (v - nint(v))**2
        besselk        = 3.747025974440712 + 6.166027005560792 * (v - nint(v))

        do j = 1, (nint(v) - 1)
            besselk_plus1 = besselk_minus1 + 8.0*skewness*besselk

            skewness = skewness + 1.0
            besselk_minus1 = besselk
            besselk = besselk_plus1
        end do
    end if
end function Kv


!------------------------------------------------------------------------------
! function Iv(v)
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
elemental function Iv(v) result(f)
    real(8) :: f
    real(8), intent(in) :: v

    real(8) :: term
    integer :: m

    term = 1.0 / (8.0**v * gamma(v + 1.0))
    f = term

    do m = 1, 6
        term = term / (64.0 * m * (v+m))
        f = f + term
    end do
end function Iv


!==============================================================================
end module shash_module
!==============================================================================
