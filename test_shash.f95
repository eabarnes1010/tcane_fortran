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
subroutine test_shash()
    write(*,*) '=========='
    write(*,*) 'test_shash'
    write(*,*) '=========='

    call test_shash_pdf_scalar()
    call test_shash_pdf_elemental()
    call test_shash_pdf_x_array()

    call test_shash_cdf_scalar()
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

    call test_shash_summary_statistics()
end subroutine test_shash


!------------------------------------------------------------------------------
subroutine test_shash_pdf_scalar()
    use shash_module
    implicit none

    real(8) :: x, loc, sigma, skewness, tailweight
    real(8) :: computed, truth

    x          = 4.5
    loc        = 1.5
    sigma      = 2.0
    skewness   = 1.5
    tailweight = 1.5

    computed = pdf(x, loc, sigma, skewness, tailweight)
    truth    = 0.08460783

    write(*,*) 'test_shash_pdf_scalar:           max_abs_error = ', abs(computed-truth)
end subroutine test_shash_pdf_scalar


!------------------------------------------------------------------------------
subroutine test_shash_pdf_elemental()
    use shash_module
    implicit none

    real(8), dimension(4) :: x, loc, sigma, skewness, tailweight
    real(8), dimension( size(x) ) :: computed, truth

    x          = [0.5, 1.0, 1.5, 2.0]
    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = pdf(x, loc, sigma, skewness, tailweight)
    truth    = [0.48394145, 0.35682481, 0.21718473, 0.1577138]

    write(*,*) 'test_shash_pdf_elemental:        max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_pdf_elemental


!------------------------------------------------------------------------------
subroutine test_shash_pdf_x_array()
    use shash_module
    implicit none

    real(8) :: x(4), loc, sigma, skewness, tailweight
    real(8), dimension( size(x) ) :: computed, truth

    x          = [1.5, 2.5, 3.5, 4.5]
    loc        = 1.5
    sigma      = 2.0
    skewness   = 1.5
    tailweight = 1.5

    computed = pdf(x, loc, sigma, skewness, tailweight)
    truth    = [0.069731, 0.16752543, 0.11880247, 0.08460783]

    write(*,*) 'test_shash_pdf_x_array:          max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_pdf_x_array


!------------------------------------------------------------------------------
subroutine test_shash_cdf_scalar()
    use shash_module
    implicit none

    real(8) :: x, loc, sigma, skewness, tailweight
    real(8) :: computed, truth

    x          = 1.5
    loc        = 1.0
    sigma      = 1.5
    skewness   = 1.0
    tailweight = 0.8

    computed = cdf(x, loc, sigma, skewness, tailweight)
    truth    = 0.22218556

    write(*,*) 'test_shash_cdf_scalar:           max_abs_error = ', abs(computed-truth)
end subroutine test_shash_cdf_scalar


!------------------------------------------------------------------------------
subroutine test_shash_cdf_elemental()
    use shash_module
    implicit none

    real(8), dimension(4) :: x, loc, sigma, skewness, tailweight
    real(8), dimension( size(x) ) :: computed, truth

    x          = [0.5, 1.0, 1.5, 2.0]
    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = cdf(x, loc, sigma, skewness, tailweight)
    truth    = [0.84134475, 0.4925046, 0.22218556, 0.07596765]

    write(*,*) 'test_shash_cdf_elemental:        max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_cdf_elemental


!------------------------------------------------------------------------------
subroutine test_shash_cdf_x_array()
    use shash_module
    implicit none

    real(8), dimension(4) :: x, loc, sigma, skewness, tailweight
    real(8), dimension( size(x) ) :: computed, truth

    x          = [1.5, 2.5, 3.5, 4.5]
    loc        = 1.5
    sigma      = 2.0
    skewness   = 1.5
    tailweight = 1.5

    computed = cdf(x, loc, sigma, skewness, tailweight)
    truth    = [0.01661558, 0.1599833, 0.3035546, 0.4036644]

    write(*,*) 'test_shash_cdf_x_array:          max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_cdf_x_array


!------------------------------------------------------------------------------
subroutine test_shash_quantile()
    use shash_module
    implicit none

    real(8), dimension(4) :: pr, loc, sigma, skewness, tailweight
    real(8), dimension( size(pr) ) :: computed, truth

    pr         = [0.1, 0.3, 0.6, 0.8]
    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = quantile(pr, loc, sigma, skewness, tailweight)
    truth    = [-0.64077578, 0.49707086, 3.46847072, 15.376874]

    write(*,*) 'test_shash_quantile (x):         max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_quantile


!------------------------------------------------------------------------------
subroutine test_shash_quantile_random()
    use shash_module
    implicit none

    integer, parameter :: NTESTS = 1000
    real(8), dimension(NTESTS) :: pr, loc, sigma, skewness, tailweight
    real(8), dimension(NTESTS) :: computed, computed_cdf

    call random_number(pr)
    call random_number(loc)
    loc = 2.0*(loc - 0.5)
    call random_number(sigma)
    sigma = 2.0*sigma
    call random_number(skewness)
    skewness = 2.0*(skewness - 0.5)
    call random_number(tailweight)
    tailweight = 2.0*tailweight

    computed     = quantile(pr, loc, sigma, skewness, tailweight)
    computed_cdf = cdf(computed, loc, sigma, skewness, tailweight)

    write(*,*) 'test_shash_quantile_random (p):  max_abs_error = ', maxval(abs(computed_cdf-pr))
end subroutine test_shash_quantile_random


!------------------------------------------------------------------------------
subroutine test_shash_median()
    use shash_module
    implicit none

    real(8), dimension(4) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(loc) ) :: computed, truth

    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = median(loc, sigma, skewness, tailweight)
    truth    = [0.0, 1.0210953, 2.8640418, 5.8619227]

    write(*,*) 'test_shash_median (x):           max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_median


!------------------------------------------------------------------------------
subroutine test_shash_median_random()
    use shash_module
    implicit none

    integer, parameter :: NTESTS = 1000
    real(8), dimension(NTESTS) :: loc, sigma, skewness, tailweight
    real(8), dimension(NTESTS) :: computed, computed_cdf

    call random_number(loc)
    loc = 2.0*(loc - 0.5)
    call random_number(sigma)
    sigma = 2.0*sigma
    call random_number(skewness)
    skewness = 2.0*(skewness - 0.5)
    call random_number(tailweight)
    tailweight = 2.0*tailweight

    computed     = median(loc, sigma, skewness, tailweight)
    computed_cdf = cdf(computed, loc, sigma, skewness, tailweight)

    write(*,*) 'test_shash_median_random (p):    max_abs_error = ', maxval(abs(computed_cdf-0.5))
end subroutine test_shash_median_random


!------------------------------------------------------------------------------
subroutine test_shash_mean()
    use shash_module
    implicit none

    real(8), dimension(4) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(loc) ) :: computed, truth

    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = mean(loc, sigma, skewness, tailweight)
    truth    = [0.0, 1.2058396, 3.2700596, 9.874768]

    write(*,*) 'test_shash_mean:                 max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_mean


!------------------------------------------------------------------------------
subroutine test_shash_mode()
    use shash_module
    implicit none

    real(8), dimension(4) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(loc) ) :: computed, truth

    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = mode(loc, sigma, skewness, tailweight)
    truth    = [0.0, 0.5764171, 1.6523985, 2.304149]

    write(*,*) 'test_shash_mode:                 max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_mode


!------------------------------------------------------------------------------
subroutine test_shash_variance()
    use shash_module
    implicit none

    real(8), dimension(4) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(loc) ) :: computed, truth

    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = variance(loc, sigma, skewness, tailweight)
    truth    = [0.25, 1.3164116, 4.4789896, 105.18572]

    write(*,*) 'test_shash_variance:             max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_variance


!------------------------------------------------------------------------------
subroutine test_shash_stddev()
    use shash_module
    implicit none

    real(8), dimension(4) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(loc) ) :: computed, truth

    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = stddev(loc, sigma, skewness, tailweight)
    truth    = [0.5, 1.1473497, 2.1163623, 10.256009]

    write(*,*) 'test_shash_stddev:               max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_stddev


!------------------------------------------------------------------------------
subroutine test_shash_skew()
    use shash_module
    implicit none

    real(8), dimension(4) :: loc, sigma, skewness, tailweight
    real(8), dimension( size(loc) ) :: computed, truth

    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = skew(loc, sigma, skewness, tailweight)
    truth    = [0.0, 0.7544219, 0.7953873, 2.29442]

    write(*,*) 'test_shash_skew:                 max_abs_error = ', maxval(abs(computed-truth))
end subroutine test_shash_skew


!------------------------------------------------------------------------------
subroutine test_shash_summary_statistics()
    use shash_module
    implicit none

    real(8), dimension(4) :: loc, sigma, skewness, tailweight
    type (summary_statistics), dimension(4) :: computed

    loc        = [0.0, 0.5, 1.0, 1.5]
    sigma      = [0.5, 1.0, 1.5, 2.0]
    skewness   = [0.0, 0.5, 1.0, 1.5]
    tailweight = [1.0, 1.0, 0.8, 1.5]

    computed = compute_summary_statistics(loc, sigma, skewness, tailweight)

    write(*,*) 'test_shash_summary_statistics'
end subroutine test_shash_summary_statistics
