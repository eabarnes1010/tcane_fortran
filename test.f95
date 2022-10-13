!==============================================================================
! test.f95                                                     (12.10.22)
!
! Test the FORTRAN-based distribution utility functions.
!
! Notes
! -----
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
!==============================================================================
PROGRAM main
    WRITE(*,*)
    CALL test_shash()
    WRITE(*,*)
    CALL test_bivariate_normal()
END PROGRAM main

