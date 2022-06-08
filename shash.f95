!==============================================================================
! shash.f95                                                          (06.06.22)
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
! Notes
! -----
! * The sinh-arcsinh normal distribution was defined in [1]. A more accessible,
! though less comprehensive, presentation is given in [2].
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
PRIVATE

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

!---------------------------------------------------------------------------
! Mathematical constants
!---------------------------------------------------------------------------
REAL(8), PARAMETER :: ONE_OVER_SQRT_TWO    = 0.7071067811865475244008444
REAL(8), PARAMETER :: ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399461

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
! FUNCTION median(mu, sigma, nu, tau)
!
! Compute the distribution median.
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
! https://www.johndcook.com/blog/normal_cdf_inverse/.
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
! FUNCTION jones_pewsey_pq(q)
!
! Compute the Jones and Pewsey factor Pq as defined on page 764 of Jones and
! Pewsey (2009).
!
! Arguments
! ---------
! q : real scalar.
!
! Returns
! -------
! p : Jones and Prewsey Pq factor.
!
!
! Notes
! -----
! * The Pq factor is used in the computation of the moments for a
! sinh-arcsinh normal (SHASH) random variable.
!
! * The Jones and Pewsey Pq factor is defined on page 764 of [1].
!
!   $P_q = \frac{e^{1/4}}{(8 \pi)^{1/2}} \left[ K_{(q+1)/2}(1/4) + K_{(q-1)/2}(1/4) \right]$
!
! where $K_{v}(1/4)$ is the modified Bessel function of the second kind.
!
! * Unfortunately, the modified Bessel function of the second kind is not
! available in the standard Fortran library. As such, we much carryout these
! computations.
!
! * We compute v = (q-1)/2, and note that (q+1)/2 = v + 1. To compute Pq we
! must compute $K_{v}(1/4)$ and $K_{v+1}(1/4)$.
!
! * Section I.2 of [2], and Section 6.5 of [3], inspired our algorithm for the
! computation of the $K_{v}(1/4)$ and $K_{v+1}(1/4)$. However, our problem is
! much simpler, since the argument is fixed at z = 1/4.
!
! We compute vo where v = vo + n, n is an integer and vo is limited
! to -0.5 < vo < 0.5. We then compute $K_{vo}(1/4)$ and $K_{vo+1}(1/4)$
! using a 15th-order Chebyshev polynomial approximation.  Finally, we compute
! the target $K_{v}(1/4)$ and $K_{v+1}(1/4)$ using a standard stable forward
! recursion.
!
! References
! ----------
! [1] M. C. Jones and A. Pewsey. "Sinh-arcsinh distributions", Biometrika 96.4
! October 2009, pp. 761–780. DOI: 10.1093/biomet/asp053.
!
! [2] N. Temme. "On the numerical evaluation of the modified Bessel function
! of the third kind", Journal of Computational Physics, 1975, 19, 324-337.
! DOI: 10.1016/0021-9991(75)90082-0.
!
! [3] S. Zhang and J. Jin. Computation of Special Functions. John Wiley and
! Sons, Inc., 1996. ISBN 978-0471119630.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION jones_pewsey_pq(q) RESULT(p)
    REAL(8) :: p
    REAL(8), INTENT(IN) :: q

    REAL(8) :: v, vo
    INTEGER :: n, j

    INTEGER, PARAMETER :: MAX_ORDER = 20
    REAL(8), DIMENSION(MAX_ORDER) :: K

    ! Partition v so that v = vo + n, where n is an integer
    ! and -0.5 < vo < 0.5.
    v = (q-1)/2
    vo = v
    n = 0
    DO WHILE (vo > 0.5)
        vo = vo - 1.0
        n = n + 1
    END DO

    ! Initialize the forward recursion on v using the 12th-order
    ! polynomial approximation.
    K(1) = limited_besselk(vo);
    K(2) = limited_besselk(vo + 1.0);

    ! Compute K_{v}(1/4) and K_{v+1}(1/4) using the standard forward
    ! recursion formula for the order. See Equation (1.9) of [2] or
    ! Equation (6.1.23) of [3].
    v = vo + 1.0
    DO j = 3, n+2
        K(j) = 8.0*v*K(j-1) + K(j-2)
        v = v + 1.0
    END DO

    ! Compute the Jones and Pewsey factor Pq. See page 764 of Jones
    ! and Pewsey (2009). The strange constant in front, 0.25612..., is
    ! $\frac{e^{1/4}}{(8\pi)^{1/2}}$ computed using a high-precision
    ! calculator.
    p = 0.25612601391340369863537463_8 * (K(n+2) + K(n+1))
END FUNCTION jones_pewsey_pq


!------------------------------------------------------------------------------
! FUNCTION limited_besselk(v)
!
! Compute $K_{v}(1/4)$ for -0.5 < v < 1.5 using a 15th-order Chebyshev
! polynomial approximation.
!
! Arguments
! ---------
! v :: real scalar in the range -0.5 < v < 1.5.
!   The range-limited order at which to evaluate $K_{v}(1/4)$.
!
! Returns
! -------
! f : real scalar
!   The approximate value of $K_{v}(1/4)$.
!
! Notes
! -----
! * This function uses 15th-order Chebyshev polynomial minimax fit for
! $K_{v}(1/4)$ over the range -0.5 < v < 1.5.
!
! * The polynomial coefficients were precomputed to minimize the maximum
! absolute error over the range -0.5 < v < 1.5. The resulting maximum
! absolute error error is 6.68e-11.
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION limited_besselk(v) RESULT(f)
    REAL(8) :: f
    REAL(8), INTENT(IN) :: v

    REAL(8) :: x, D1, D2, S
    INTEGER :: j

    ! Precomputed coefficients for a 15th-order Chebyshev polynomial minimax
    ! fit of $K_{v}(1/4)$ over the range -0.5 < v < 1.5.
    REAL(8), PARAMETER :: C(16) = (/&
        3.703498483049486e+00_8, &
        3.291364878445456e+00_8, &
        1.940643036786937e+00_8, &
        5.662030547496981e-01_8, &
        2.003994693199613e-01_8, &
        4.459824835968565e-02_8, &
        1.150959269394563e-02_8, &
        2.093789904913530e-03_8, &
        4.299445769978876e-04_8, &
        6.657921299844640e-05_8, &
        1.143004469426901e-05_8, &
        1.547568909556892e-06_8, &
        2.293008750401813e-07_8, &
        2.765961264284608e-08_8, &
        3.614997877748472e-09_8, &
        4.019025718452231e-10_8 /)

    ! Rescale the interval from -0.5 : 1.5 to -1 : 1.
    x = v - 0.5

    ! Evaluate the 15th-order Chebyshev polynomial using Clenshaw recurrence.
    DO j = 16, 2, -1
        S = D1;
        D1 = 2*x*D1 - D2 + c(j);
        D2 = S;
    END DO
    f = x*D1 - D2 + c(1);
END FUNCTION limited_besselk

!==============================================================================
END MODULE shash_module
!==============================================================================
