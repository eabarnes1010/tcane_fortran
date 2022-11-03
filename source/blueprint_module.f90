!==============================================================================
! MODULE: blueprints_module
!
! Define the components of the artificial neural network to be passed from the
! Python-based training code (Tensorflow) to the stand-alone, Fortran-based,
! evaluation code.
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
! * 03 November 2022
!
!==============================================================================
module blueprint_module
   use read_json_module

   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private

   public :: read_blueprint

   !-----------------------------------
   ! Module parameters
   !-----------------------------------
   integer, parameter :: EUNIT = 11

   !-----------------------------------
   ! TYPE: InputLayerTraits
   !-----------------------------------
   type, public :: InputLayerTraits
      real(8), dimension(:), allocatable :: mean
      real(8), dimension(:), allocatable :: std
   end type InputLayerTraits

   !-----------------------------------
   ! TYPE: HiddenLayerTraits
   !-----------------------------------
   type, public :: HiddenLayerTraits
      real(8), dimension(:,:), allocatable :: weights
      real(8), dimension(:),   allocatable :: bias
      character(len=16)                    :: activation
   end type HiddenLayerTraits

   !-----------------------------------
   ! TYPE: OutputChannelTraits
   !-----------------------------------
   type, public :: OutputChannelTraits
      real(8), dimension(:), allocatable :: weights
      real(8)                            :: bias
      character(len=16)                  :: transformation
      real(8)                            :: mean
      real(8)                            :: std
   end type OutputChannelTraits

   !-----------------------------------
   ! TYPE: TestCases
   !-----------------------------------
   type, public :: TestCases
      integer                              :: n_test
      real(8), dimension(:,:), allocatable :: x_test
      real(8), dimension(:,:), allocatable :: predictions
   end type TestCases


   !-----------------------------------
   ! TYPE: Blueprint
   !-----------------------------------
   type, public :: Blueprint
      integer                                              :: n_input
      integer                                              :: n_hidden
      integer                                              :: n_output

      character(len=16),         dimension(:), allocatable :: input_names
      character(len=16),         dimension(:), allocatable :: output_names

      type(InputLayerTraits)                               :: input_traits
      type(HiddenLayerTraits),   dimension(:), allocatable :: hidden_traits
      type(OutputChannelTraits), dimension(:), allocatable :: output_traits

      type(TestCases),                         allocatable :: test_cases
   end type Blueprint

contains

   !-----------------------------------
   ! FUNCTION: read_blueprint
   !-----------------------------------
   function read_blueprint(filename) result(output)
      character(*), intent(in) :: filename
      type(Blueprint), allocatable :: output

      integer alloc_stat, io_stat, i

      allocate(output, stat=alloc_stat)
      if (alloc_stat /= 0) error stop 'allocation error in <read_blueprint #1>'

      open(unit=EUNIT, file=filename, status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) error stop 'file open error in <read_blueprint>'

      ! Skip the { on the first line.
      call skip_line(eunit)

      ! Read the general traits.
      output%n_input  = read_integer(eunit)
      output%n_hidden = read_integer(eunit)
      output%n_output = read_integer(eunit)

      output%input_names  = read_list(eunit, output%n_input)
      output%output_names = read_list(eunit, output%n_output)

      ! Read the input traits.
      call skip_line(eunit)
      output%input_traits = read_input_traits(eunit)
      call skip_line(eunit)

      ! Read the hidden traits one layer at a time.
      allocate(output%hidden_traits(output%n_hidden), stat=alloc_stat)
      if (alloc_stat /= 0) error stop 'allocation error in <read_blueprint #2>'

      call skip_line(eunit)
      do i = 1, output%n_hidden
         call skip_line(eunit)
         output%hidden_traits(i) = read_hidden_traits(eunit)
         call skip_line(eunit)
      end do
      call skip_line(eunit)

      ! Read the output traits one channel at a time.
      allocate(output%output_traits(output%n_output), stat=alloc_stat)
      if (alloc_stat /= 0) error stop 'allocation error in <read_blueprint #3>'

      call skip_line(eunit)
      do i = 1, output%n_output
         call skip_line(eunit)
         output%output_traits(i) = read_output_traits(eunit)
         call skip_line(eunit)
      end do
      call skip_line(eunit)

      ! Read the test cases.
      call skip_line(eunit)
      output%test_cases = read_test_cases(eunit, output%n_input, output%n_output)
      call skip_line(eunit)

      close(unit=EUNIT)
   end function read_blueprint

   !-----------------------------------
   ! FUNCTION: read_input_traits
   !-----------------------------------
   function read_input_traits(eunit) result(output)
      integer :: eunit
      type(InputLayerTraits), allocatable :: output

      integer alloc_stat, n_input

      allocate(output, stat=alloc_stat)
      if (alloc_stat /= 0) error stop 'allocation error in <read_input_traits>'

      n_input     = read_integer(eunit)
      output%mean = read_vector(eunit, n_input)
      output%std  = read_vector(eunit, n_input)
   end function read_input_traits

   !-----------------------------------
   ! FUNCTION: read_hidden_traits
   !-----------------------------------
   function read_hidden_traits(eunit) result(output)
      integer :: eunit
      type(HiddenLayerTraits), allocatable :: output

      integer alloc_stat, n_input, n_output

      allocate(output, stat=alloc_stat)
      if (alloc_stat /= 0) error stop 'allocation error in <read_hidden_traits>'

      n_input  = read_integer(eunit)
      n_output = read_integer(eunit)

      output%weights    = read_matrix(eunit, n_input, n_output)
      output%bias       = read_vector(eunit, n_output)
      output%activation = read_string(eunit)
   end function read_hidden_traits

   !-----------------------------------
   ! FUNCTION: read_output_traits
   !-----------------------------------
   function read_output_traits(eunit) result(output)
      integer :: eunit
      type(OutputChannelTraits), allocatable :: output

      integer alloc_stat, n_input

      allocate(output, stat=alloc_stat)
      if (alloc_stat /= 0) error stop 'allocation error in <read_output_traits>'

      n_input  = read_integer(eunit)

      output%weights        = read_vector(eunit, n_input)
      output%bias           = read_real(eunit)
      output%transformation = read_string(eunit)
      output%mean           = read_real(eunit)
      output%std            = read_real(eunit)
   end function read_output_traits

   !-----------------------------------
   ! FUNCTION: read_test_cases
   !-----------------------------------
   function read_test_cases(eunit, n_input, n_output) result(output)
      integer :: eunit, n_input, n_output
      type(TestCases), allocatable :: output

      integer alloc_stat

      allocate(output, stat=alloc_stat)
      if (alloc_stat /= 0) error stop 'allocation error in <read_test_cases>'

      output%n_test      = read_integer(eunit)
      output%x_test      = read_matrix(eunit, output%n_test, n_input)
      output%predictions = read_matrix(eunit, output%n_test, n_output)
   end function read_test_cases

end module blueprint_module
