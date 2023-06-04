!==============================================================================
! MODULE: input_layer_module
!
! Define and deploy the input layer with normalization.
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
! * 25 October 2022
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

      if (size(mean) /= size(std)) then
         print *, size(mean), size(std)
         error stop "mismatched mean and std arrays in <initialize_InputLayer>"
      end if

      if (any(std <= 0.0_8)) then
         error stop "invalid std array in <initialize_InputLayer>"
      end if

      this%n_input = size(mean)

      if (allocated(this%mean)) deallocate(this%mean)
      
      !apparently some gfortran compilers do not support the
      !source= specifier in the allocate statement, even with using -std=f2003 or f2008.
      !Thus replace:
      !allocate(this%mean, source=mean, stat=alloc_stat)
      !by: (similar changes in other places that use source= specifier)
      allocate(this%mean(this%n_input), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "allocation error in <in <initialize_InputLayer #1>"
      this%mean = mean

      if (allocated(this%std)) deallocate(this%std)
      !allocate(this%std, source=std, stat=alloc_stat)
      allocate(this%std(this%n_input), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "allocation error in <in <initialize_InputLayer #2>"
      this%std = std
   end subroutine initialize_InputLayer

   !-----------------------------------
   ! FUNCTION: employ_InputLayer
   !-----------------------------------
   function employ_InputLayer(this, input) result(output)
      type(InputLayer),                intent(in) :: this
      real(8), dimension(:),           intent(in) :: input
      real(8), dimension(this%n_input)            :: output

      if (size(input) /= this%n_input)  error stop "incorrect input array size in <employ_InputLayer>"

      output = (input - this%mean)/this%std
   end function employ_InputLayer

end module input_layer_module
