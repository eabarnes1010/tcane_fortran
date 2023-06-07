!==============================================================================
! MODULE: output_channel_module
!
! Define and deploy the output channels with weights, bias, activation,
! transformation, and scaling.
!
! Notes
! -----
! * An output channel has a single neuron.  This neuron is fully-connected to
!   the last hidden layer. Thus, the output channel has a vector of weights,
!   a scalar bias, and an activation function.
!
!   The scalar output from the neuron is then fed into a transformation
!   function. For channels that must yield strictly positive values (e.g. a
!   standard deviation), this transformation will be a Softmax. For channels
!   that are bounded between plus/minus one (e.g. a correlation), this trans-
!   formation will be a Tanh. For unconstrained cases, this transformation
!   will be a simple pass-through using a Linear.
!
!   The scalar output from the transformation is then fed into rescaling.  The
!   rescaling is the opposite of normalization. The rescaling multiplies by the
!   std and adds the mean.
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

      if (.not. (is_valid_transformation(transformation))) then
         print *, transformation
         error stop "invalid transformation in <initialize_OutputChannel>"
      end if

      if (std < 0.0_8) then
         print *, std
         error stop "invalid std in <initialize_OutputChannel>"
      end if

      this%n_input  = size(weights)

      if (allocated(this%weights)) deallocate(this%weights)
      !allocate(this%weights, source=weights, stat=alloc_stat)
      allocate(this%weights(this%n_input), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "allocation error in <in <initialize_OutputChannel>"
      this%weights = weights

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

      if (size(input) /= this%n_input)  error stop "incorrect input array size in <employ_OutputChannel>"

      throughput = dot_product(this%weights, input) + this%bias
      throughput = apply_transformation(this%transformation, throughput)
      output     = this%std*throughput + this%mean
   end function employ_OutputChannel

end module output_channel_module
