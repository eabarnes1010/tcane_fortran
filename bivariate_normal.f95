!==============================================================================
! bivariate_normal.f95                                               (10.12.22)
!
! FORTRAN-based bivariate normal distribution utility functions.
!
! Public functions
! ----------------
! pdf(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
!   Compute the bivariate normal probability density function (pdf).
!   * u, v, mu_u, mu_v, sigma_u, sigma_v, and rho may be scalars.
!   * u, v, mu_u, mu_v, sigma_u, sigma_v, and rho may be commensurate
!       arrays.
!   * u and v may be commensurate arrays, with scalar mu_u, mu_v,
!       sigma_u, sigma_v, and rho.
!
! cdf(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
!   Compute the Mahalanobis cumulative distribution function (cdf).
!   * u, v, mu_u, mu_v, sigma_u, sigma_v, and rho may be scalars.
!   * u, v, mu_u, mu_v, sigma_u, sigma_v, and rho may be commensurate
!       arrays.
!   * u and v may be commensurate arrays, with scalar mu_u, mu_v,
!       sigma_u, sigma_v, and rho
!
! inv_cdf(pr, mu_u, mu_v, sigma_u, sigma_v, rho, u, v)
!   Compute arrays of u and v coordinates defining the Mahalanobis ellipse
!   that captures a probability of pr.
!
! Notes
! -----
! * This code was tested using:
!   GNU Fortran (MinGW-W64 x86_64-ucrt-posix-seh, built by Brecht Sanders) 12.2.0
!
! References
! ----------
! [1] Michael Bensimhoun. (2009, June). N-dimension Cumulative Function
!     and Other Useful Facts About Gaussians and Normal Densities.
!     https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf
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
MODULE bivariate_normal_module
IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
! Interface block
!------------------------------------------------------------------------------
INTERFACE pdf
    MODULE PROCEDURE pdf_elemental
    MODULE PROCEDURE pdf_uv_array
END INTERFACE

INTERFACE cdf
    MODULE PROCEDURE mahalanobis_cdf_elemental
    MODULE PROCEDURE mahalanobis_cdf_uv_array
END INTERFACE

INTERFACE invcdf
    MODULE PROCEDURE mahalanobis_inv_cdf
END INTERFACE


!------------------------------------------------------------------------------
! Visibility
!------------------------------------------------------------------------------
PUBLIC pdf
PUBLIC cdf
PUBLIC invcdf

!---------------------------------------------------------------------------
! Mathematical constants
!---------------------------------------------------------------------------
REAL(8), PARAMETER :: PI                   = 3.1415926535897932384626434_8
REAL(8), PARAMETER :: TWO_PI               = 6.2831853071795864769252868_8
REAL(8), PARAMETER :: PI_OVER_TWO          = 1.5707963267948966192313217_8
REAL(8), PARAMETER :: ONE_OVER_SQRT_TWO    = 0.7071067811865475244008444_8
REAL(8), PARAMETER :: ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399461_8

CONTAINS

!------------------------------------------------------------------------------
! FUNCTION pdf_elemental(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
!
! Compute the bivariate normal probability density function (pdf).
!
! Arguments
! ---------
! u : real scalar
!   u value for evaluation of pdf.
!
! v : real scalar
!   v value for evaluation of pdf.
!
! mu_u : real scalar
!   mean of u.
!
! mu_v : real scalar
!   mean of v.
!
! sigma_u : real scalar, sigma_u > 0.
!   standard deviation of u.
!
! sigma_v : real scalar, sigma_v > 0.
!   standard deviation of v.
!
! rho : real scalar, -1 < rho < 1.
!   correlation between u and v.
!
! Returns
! -------
! pdf : real scalar
!   The computed bivariate normal probability density function
!   evaluated at (u, v).
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
! * The equation for the returned pdf comes from the the bottom of
!   Page 2 of [1].
!
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION pdf_elemental(u, v, mu_u, mu_v, sigma_u, sigma_v, rho) RESULT(pdf)
    REAL(8), INTENT(IN) :: u, v, mu_u, mu_v, sigma_u, sigma_v, rho
    REAL(8) :: pdf

    REAL(8) :: s, t, c, r_sqr

    s = (u - mu_u)/sigma_u
    t = (v - mu_v)/sigma_v
    c = 1.0/(2.0*PI * sigma_u*sigma_v * SQRT(1.0 - rho**2))
    r_sqr = 1.0 / (1.0 - rho**2) * (s*s - 2.0*rho*s*t + t*t)

    pdf = c * EXP(-r_sqr/2.0)
END FUNCTION pdf_elemental


!------------------------------------------------------------------------------
! FUNCTION pdf_uv_array(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
!
! Compute the bivariate normal probability density function (pdf).
!
! Arguments
! ---------
! u : real array, size(u) = size(v)
!   u value for evaluation of pdf.
!
! v : real scalar, size(u) = size(v)
!   v value for evaluation of pdf.
!
! mu_u : real scalar
!   mean of u.
!
! mu_v : real scalar
!   mean of v.
!
! sigma_u : real scalar, sigma_u > 0.
!   standard deviation of u.
!
! sigma_v : real scalar, sigma_v > 0.
!   standard deviation of v.
!
! rho : real scalar, -1 < rho < 1.
!   correlation between u and v.
!
! Returns
! -------
! pdf : real array
!   The computed bivariate normal probability density function
!   evaluated at (u, v). pdf is the same size as u and v.
!
! Notes
! -----
! * The equation for the returned pdf comes from the the bottom of
!   Page 2 of [1].
!
!------------------------------------------------------------------------------
PURE FUNCTION pdf_uv_array(u, v, mu_u, mu_v, sigma_u, sigma_v, rho) RESULT(pdf)
    REAL(8), INTENT(IN) :: u(:), v(:)
    REAL(8), INTENT(IN) :: mu_u, mu_v, sigma_u, sigma_v, rho
    REAL(8), DIMENSION( SIZE(u) ) :: pdf

    REAL(8), DIMENSION( SIZE(u)) :: s, t, c, r_sqr

    s = (u - mu_u)/sigma_u
    t = (v - mu_v)/sigma_v
    c = 1.0/(2.0*PI * sigma_u * sigma_v * SQRT(1.0 - rho**2))
    r_sqr = 1.0 / (1.0 - rho**2) * (s*s - 2.0*rho*s*t + t*t)

    pdf = c * EXP(-r_sqr/2.0)
END FUNCTION pdf_uv_array


!------------------------------------------------------------------------------
! FUNCTION mahalanobis_cdf_elemental(u, v, mu_u, mu_v, sigma_u, sigma_v, rho)
!
! Compute the bivariate normal Mahalanobis cumulative distribution
! function (cdf).  That is, the probability that a point falls inside the
! ellipse defined by the Mahalanobis distance given by (u, v).
!
! Arguments
! ---------
! u : real scalar
!   u value for evaluation of pdf.
!
! v : real scalar
!   v value for evaluation of pdf.
!
! mu_u : real scalar
!   mean of u.
!
! mu_v : real scalar
!   mean of v.
!
! sigma_u : real scalar, sigma_u > 0.
!   standard deviation of u.
!
! sigma_v : real scalar, sigma_v > 0.
!   standard deviation of v.
!
! rho : real scalar, -1 < rho < 1.
!   correlation between u and v.
!
! Returns
! -------
! cdf : real scalar
!   The computed bivariate normal Mahalanobis cumulative distribution
!   function evaluated at (u, v).
!
! Notes
! -----
! * This function is ELEMENTAL, so it may be called with scalar arguments. It
!   may also be called with array arguments as long as all of the arrays are
!   commensurate.
!
! * The equations for the Mahalanobis distance, r_sqr, comes from
!   the bottom of Page 2 in [1].
!
! * The equation for the returned cdf comes from the the middle of
!   Page 4 of [1].
!
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION mahalanobis_cdf_elemental(u, v, mu_u, mu_v, sigma_u, sigma_v, rho) RESULT(cdf)
    REAL(8), INTENT(IN) :: u, v, mu_u, mu_v, sigma_u, sigma_v, rho
    REAL(8) :: cdf

    REAL(8) :: s, t, r_sqr

    s = (u - mu_u)/sigma_u
    t = (v - mu_v)/sigma_v
    r_sqr = 1.0 / (1.0 - rho**2) * (s*s - 2.0*rho*s*t + t*t)

    cdf = 1.0 - EXP(-r_sqr/2.0)
END FUNCTION mahalanobis_cdf_elemental


!------------------------------------------------------------------------------
! FUNCTION cdf_uv_array(u, v, mu_u, mu_v, sigma_u, sigma_v, rho, r, x, y)
!
! Compute the bivariate normal Mahalanobis cumulative distribution
! function (cdf).  That is, the probability that a point falls inside the
! ellipse defined by the Mahalanobis distance given by (u, v).
!
! Arguments
! ---------
! u : real array, size(u) = size(v)
!   u value for evaluation of pdf.
!
! v : real scalar, size(u) = size(v)
!   v value for evaluation of pdf.
!
! mu_u : real scalar
!   mean of u.
!
! mu_v : real scalar
!   mean of v.
!
! sigma_u : real scalar, sigma_u > 0.
!   standard deviation of u.
!
! sigma_v : real scalar, sigma_v > 0.
!   standard deviation of v.
!
! rho : real scalar, -1 < rho < 1.
!   correlation between u and v.
!
! Returns
! -------
! pdf : real array
!   The computed bivariate normal Mahalanobis cumulative distribution
!   function evaluated at (u, v).
!
! Notes
! -----
! * The equations for the Mahalanobis distance, r_sqr, comes from
!   the bottom of Page 2 in [1].
!
! * The equation for the returned cdf comes from the the middle of
!   Page 4 of [1].
!
!------------------------------------------------------------------------------
PURE FUNCTION mahalanobis_cdf_uv_array(u, v, mu_u, mu_v, sigma_u, sigma_v, rho) RESULT(cdf)
    REAL(8), INTENT(IN) :: u(:), v(:)
    REAL(8), INTENT(IN) :: mu_u, mu_v, sigma_u, sigma_v, rho
    REAL(8), DIMENSION( SIZE(u) ) :: cdf

    REAL(8), DIMENSION( SIZE(u)) :: s, t, r_sqr

    s = (u - mu_u)/sigma_u
    t = (v - mu_v)/sigma_v
    r_sqr = 1.0 / (1.0 - rho**2) * (s*s - 2.0*rho*s*t + t*t)

    cdf = 1.0 - EXP(-r_sqr/2.0)
END FUNCTION mahalanobis_cdf_uv_array


!------------------------------------------------------------------------------
! SUBROUTINE mahalanobis_inv_cdf(pr, mu_u, mu_v, sigma_u, sigma_v, rho, u, v)
!
! Compute arrays of u and v coordinates defining the Mahalanobis ellipse that
! captures a probability of pr.
!
! Arguments
! ---------
! pr : real scalar
!   probability of capture.
!
! mu_u : real scalar
!   mean of u.
!
! mu_v : real scalar
!   mean of v.
!
! sigma_u : real scalar, sigma_u > 0.
!   standard deviation of u.
!
! sigma_v : real scalar, sigma_v > 0.
!   standard deviation of v.
!
! rho : real scalar, -1 < rho < 1.
!   correlation between u and v.
!
! u : real array
!   returned array of u-coordinates defining the ellipse
!
! v : real array
!   returned array of u-coordinates defining the ellipse
!
! Notes
! -----
! * The equation for r comes the bottom of Page 4 of [1].  The equations for
!   u and v comes from the bottom of Page 2 of [1].
!
! * The array of returned u and v coordinates are given in counter-clockwise
!   order and the first and last points coincide.
!
!------------------------------------------------------------------------------
SUBROUTINE mahalanobis_inv_cdf(pr, mu_u, mu_v, sigma_u, sigma_v, rho, u, v)
    REAL(8), INTENT(IN) :: pr, mu_u, mu_v, sigma_u, sigma_v, rho

    INTEGER :: I
    INTEGER, PARAMETER :: NPOINTS = 1001
    REAL(8), DIMENSION(NPOINTS), PARAMETER :: THETA = (/ (REAL(I)*TWO_PI/(NPOINTS-1), I = 0, NPOINTS-1) /)
    REAL(8), DIMENSION(NPOINTS), INTENT(OUT) :: u, v

    REAL(8) :: r

    r = SQRT(-2.0 * (LOG(1.0 - pr)))
    u = r * sigma_u * COS(THETA) + mu_u
    v = r * sigma_v * (rho * COS(THETA) + SQRT(1.0 - rho**2) * SIN(THETA)) + mu_v
END SUBROUTINE mahalanobis_inv_cdf


!==============================================================================
END MODULE bivariate_normal_module
!==============================================================================
