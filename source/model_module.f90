!==============================================================================
! MODULE: model_module
!
! Define the model with an input layer, multiple hidden layers, and an output
! layer.
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
module model_module
   use blueprint_module
   use hidden_layer_module
   use input_layer_module
   use output_channel_module

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: initialize
   public :: employ

   !-----------------------------------
   ! TYPE: Model
   !-----------------------------------
   type, public :: Model
      private

      integer :: n_input
      integer :: n_hidden
      integer :: n_output

      character(len=16),   dimension(:), allocatable :: input_names
      character(len=16),   dimension(:), allocatable :: output_names

      type(InputLayer)                               :: input_layer
      type(HiddenLayer),   dimension(:), allocatable :: hidden_layers
      type(OutputChannel), dimension(:), allocatable :: output_channels
   end type Model

   !-----------------------------------
   ! Interfaces
   !-----------------------------------
   interface initialize
      module procedure initialize_Model
   end interface initialize

   interface employ
      module procedure employ_Model
   end interface

contains

   !-----------------------------------
   ! SUBROUTINE: initialize_Model
   !-----------------------------------
   subroutine initialize_Model(this, plans)
      type(Model),     intent(inout) :: this
      type(Blueprint), intent(in)    :: plans

      integer :: i, alloc_stat

      this%n_input  = plans%n_input
      this%n_hidden = plans%n_hidden
      this%n_output = plans%n_output

      this%input_names  = plans%input_names
      this%output_names = plans%output_names

      call initialize( &
         this%input_layer, &
         plans%input_traits%mean, &
         plans%input_traits%std &
      )

      if (allocated(this%hidden_layers)) deallocate(this%hidden_layers)
      allocate(this%hidden_layers(this%n_hidden), stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation error in <in <initialize_Model #1>"

      do i = 1, this%n_hidden
         call initialize( &
            this%hidden_layers(i), &
            plans%hidden_traits(i)%weights, &
            plans%hidden_traits(i)%bias, &
            plans%hidden_traits(i)%activation &
         )
      end do

      if (allocated(this%output_channels)) deallocate(this%output_channels)
      allocate(this%output_channels(this%n_output), stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation error in <in <initialize_Model #2>"

      do i = 1, this%n_output
         call initialize( &
            this%output_channels(i), &
            plans%output_traits(i)%weights, &
            plans%output_traits(i)%bias, &
            plans%output_traits(i)%transformation, &
            plans%output_traits(i)%mean, &
            plans%output_traits(i)%std &
         )
      end do
   end subroutine initialize_Model

   !-----------------------------------
   ! FUNCTION: employ_Model
   !-----------------------------------
   function employ_Model(this, input) result(output)
      type(Model),                      intent(in) :: this
      real(8), dimension(this%n_input), intent(in) :: input
      real(8), dimension(this%n_output)            :: output

      integer :: i
      real(8), dimension(:), allocatable :: throughput
      allocate(throughput(0))

      if (size(input) /= this%n_input) stop "incorrect input array size in <employ_Model>"

      throughput = employ(this%input_layer, input)

      do i = 1, this%n_hidden
         throughput = employ(this%hidden_layers(i), throughput)
      end do

      do i = 1, this%n_output
         output(i) = employ(this%output_channels(i), throughput)
      end do
   end function employ_Model

end module model_module
