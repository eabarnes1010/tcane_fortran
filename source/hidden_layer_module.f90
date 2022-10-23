!==============================================================================
! MODULE: hidden_layer_module
!
! Define the dense, hidden, layers with neurons.
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
! * 24 October 2022
!
!==============================================================================
module hidden_layer_module
   use activation_module

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: initialize
   public :: employ

   !-----------------------------------
   ! TYPE: HiddenLayer
   !-----------------------------------
   type, public :: HiddenLayer
      private

      integer                                        :: n_input
      integer                                        :: n_output

      real(8),          dimension(:, :), allocatable :: weights
      real(8),          dimension(:),    allocatable :: bias
      character(len=12)                              :: activation
   end type HiddenLayer

   !-----------------------------------
   ! Interfaces
   !-----------------------------------
   interface initialize
      module procedure initialize_HiddenLayer
   end interface initialize

   interface employ
      module procedure employ_HiddenLayer
   end interface

contains

   !-----------------------------------
   ! SUBROUTINE: initialize_HiddenLayer
   !-----------------------------------
   subroutine initialize_HiddenLayer(this, weights, bias, activation)
      type(HiddenLayer),         intent(inout) :: this
      real(8), dimension(:, :), intent(in)    :: weights
      real(8), dimension(:),    intent(in)    :: bias
      character(12),            intent(in)    :: activation

      integer :: alloc_stat

      if (size(weights, 1) /= size(bias)) then
         print *, size(weights)
         print *, size(bias)
         stop "incompatible weights and bias arrays in <initialize_HiddenLayer>"
      end if

      this%n_input  = size(weights, 2)
      this%n_output = size(bias)

      if (allocated(this%weights)) deallocate(this%weights)
      allocate(this%weights, source=weights, stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation error in <in <initialize_HiddenLayer #1>"

      if (allocated(this%bias)) deallocate(this%bias)
      allocate(this%bias, source=bias, stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation error in <in <initialize_HiddenLayer #2>"

      this%activation = activation
   end subroutine initialize_HiddenLayer

   !-----------------------------------
   ! FUNCTION: employ_HiddenLayer
   !-----------------------------------
   function employ_HiddenLayer(this, input) result(output)
      type(HiddenLayer),                         intent(in) :: this
      real(8),          dimension(:),            intent(in) :: input
      real(8),          dimension(this%n_output)            :: output

      real(8), dimension(this%n_output) :: throughput

      if (size(input) /= this%n_input)  stop "incorrect input array size in <employ_HiddenLayer>"

      throughput = matmul(transpose(this%weights), input) + this%bias
      output     = apply_activation(this%activation, throughput)
   end function employ_HiddenLayer

end module hidden_layer_module
