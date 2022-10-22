!==============================================================================
! MODULE: output_channel_module
!
! Define the output channel with weights, bias, activation, transformation,
! and scaling.
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
! * 22 October 2022
!
!==============================================================================
module output_channel_module
   use activation_module
   use transformation_module

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: initialize
   public :: employ

   !-----------------------------------
   ! TYPE: OutputChannel
   !-----------------------------------
   type, public :: OutputChannel
      private

      integer                                     :: n_input

      real(8),          dimension(:), allocatable :: weights
      real(8)                                     :: bias
      character(len=12)                           :: transformation
      real(8)                                     :: mean
      real(8)                                     :: std
   end type OutputChannel

   !-----------------------------------
   ! Interfaces
   !-----------------------------------
   interface initialize
      module procedure initialize_OutputChannel
   end interface initialize

   interface employ
      module procedure employ_OutputChannel
   end interface

contains

   !-----------------------------------
   ! SUBROUTINE: initialize_OutputChannel
   !-----------------------------------
   subroutine initialize_OutputChannel(this, weights, bias, transformation, mean, std)
      type(OutputChannel),               intent(inout) :: this
      real(8),             dimension(:), intent(in)    :: weights
      real(8),                           intent(in)    :: bias
      character(len=12),                 intent(in)    :: transformation
      real(8),                           intent(in)    :: mean
      real(8),                           intent(in)    :: std

      integer :: alloc_stat

      this%n_input  = size(weights)

      if (allocated(this%weights)) deallocate(this%weights)
      allocate(this%weights, source=weights, stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation failed"

      this%bias           = bias
      this%transformation = transformation
      this%mean           = mean
      this%std            = std
   end subroutine initialize_OutputChannel

   !-----------------------------------
   ! FUNCTION: employ_OutputChannel
   !-----------------------------------
   function employ_OutputChannel(this, input) result(output)
      type(OutputChannel),               intent(in) :: this
      real(8),             dimension(:), intent(in) :: input
      real(8)                                       :: output

      real(8) :: throughput

      if (size(input) /= this%n_input)  stop "incorrect input array size"

      throughput = dot_product(this%weights, input) + this%bias
      throughput = apply_transformation(this%transformation, throughput)
      output     = this%std*throughput + this%mean
   end function employ_OutputChannel

end module output_channel_module