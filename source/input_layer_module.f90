!==============================================================================
! MODULE: input_layer_module
!
! Define the output layer with normalization.
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
! * 21 October 2022
!
!==============================================================================
module input_layer_module
   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private
   public :: initialize
   public :: employ

   !-----------------------------------
   ! TYPE: InputLayer
   !-----------------------------------
   type, public :: InputLayer
      private

      integer                            :: n_input

      real(8), dimension(:), allocatable :: mean
      real(8), dimension(:), allocatable :: std
   end type InputLayer

   !-----------------------------------
   ! Interfaces
   !-----------------------------------
   interface initialize
      module procedure initialize_InputLayer
   end interface initialize

   interface employ
      module procedure employ_InputLayer
   end interface

contains

   !-----------------------------------
   ! SUBROUTINE: initialize_InputLayer
   !-----------------------------------
   subroutine initialize_InputLayer(this, mean, std)
      type(InputLayer),      intent(inout) :: this
      real(8), dimension(:), intent(in)    :: mean
      real(8), dimension(:), intent(in)    :: std

      integer :: alloc_stat

      if (size(mean) /= size(std)) stop "mismatched mean and std arrays"

      this%n_input = size(mean)

      if (allocated(this%mean)) deallocate(this%mean)
      allocate(this%mean, source=mean, stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation failed"

      if (allocated(this%std)) deallocate(this%std)
      allocate(this%std, source=std, stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation failed"
   end subroutine initialize_InputLayer

   !-----------------------------------
   ! FUNCTION: employ_InputLayer
   !-----------------------------------
   function employ_InputLayer(this, input) result(output)
      type(InputLayer),                intent(in) :: this
      real(8), dimension(:),           intent(in) :: input
      real(8), dimension(this%n_input)            :: output

      if (size(input) /= this%n_input)  stop "incorrect input array size"

      output = (input - this%mean)/this%std
   end function employ_InputLayer

end module input_layer_module
