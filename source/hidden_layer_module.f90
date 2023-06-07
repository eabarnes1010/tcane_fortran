!==============================================================================
! MODULE: hidden_layer_module
!
! Define and deploy the fully-connected hidden layers with weights, bias, and
! activation.
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
! Dr. Galina Chirokova
! Cooperative Institute for Research in the Atmosphere
! Colorado State University
! Fort Collins, CO, 80523
!
! Version
! -------
! * 05 June 2023
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

      if (size(weights, 2) /= size(bias)) then
         print *, size(weights, 1), size(weights, 2)
         print *, size(bias)
         error stop "incompatible weights and bias arrays in <initialize_HiddenLayer>"
      end if

      if (.not. (is_valid_activation(activation))) then
         print *, activation
         error stop "invalid activation in <initialize_HiddenLayer>"
      end if

      this%n_input    = size(weights, 1)
      this%n_output   = size(weights, 2)
      this%activation = activation

      if (allocated(this%weights)) deallocate(this%weights)
      !allocate(this%weights, source=weights, stat=alloc_stat)
      allocate(this%weights(this%n_input, this%n_output), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "allocation error in <initialize_HiddenLayer #1>"
      this%weights = weights

      if (allocated(this%bias)) deallocate(this%bias)
      !allocate(this%bias, source=bias, stat=alloc_stat)
      allocate(this%bias(this%n_output), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "allocation error in <initialize_HiddenLayer #2>"
      this%bias = bias
   end subroutine initialize_HiddenLayer

   !-----------------------------------
   ! FUNCTION: employ_HiddenLayer
   !-----------------------------------
   function employ_HiddenLayer(this, input) result(output)
      type(HiddenLayer),                         intent(in) :: this
      real(8),          dimension(:),            intent(in) :: input
      real(8),          dimension(this%n_output)            :: output

      real(8), dimension(this%n_output) :: throughput

      if (size(input) /= this%n_input) then
         print *, size(input)
         error stop "incorrect input array size in <employ_HiddenLayer>"
      end if

      throughput = matmul(input, this%weights) + this%bias
      output     = apply_activation(this%activation, throughput)
   end function employ_HiddenLayer

end module hidden_layer_module
