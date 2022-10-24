!==============================================================================
! MODULE: read_json_module
!
!
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
module read_json_module
   implicit none

   !-----------------------------------
   ! Visibility
   !-----------------------------------
   private

   public :: read_integer
   public :: read_list
   public :: read_matrix
   public :: read_real
   public :: read_string
   public :: read_vector
   public :: skip_line

   !-----------------------------------
   ! Module parameters
   !-----------------------------------
   integer,      parameter :: STRING_SIZE = 16
   integer,      parameter :: LINE_SIZE   = 1024
   character(*), parameter :: LINE_FORMAT = '(A1024)'

contains

   !-----------------------------------
   ! FUNCTION: get_line
   !-----------------------------------
   function get_line(eunit) result(line)
      integer, intent(in) :: eunit
      character(len=LINE_SIZE) :: line

      integer :: io_stat

      read(eunit, LINE_FORMAT, iostat=io_stat) line
      if (io_stat /= 0) then
         print *, '>>> ' // adjustl(trim(line))
         stop 'read error in <get_line>'
      end if

      ! print *, trim(line)
   end function get_line

   !-----------------------------------
   ! FUNCTION: read_integer
   !-----------------------------------
   function read_integer(eunit) result(output)
      integer, intent(in) :: eunit
      integer :: output

      character(len=LINE_SIZE) :: line
      integer :: io_stat, left

      line = get_line(eunit)
      left = scan(line, ':') + 1

      read(line(left:), *, iostat=io_stat) output
      if (io_stat /= 0) then
         print *, '>>> ' // adjustl(trim(line))
         stop 'read error in <get_integer>'
      end if
   end function read_integer

   !-----------------------------------
   ! FUNCTION: read_list
   !-----------------------------------
   function read_list(eunit, nrows) result(output)
      integer, intent(in) :: eunit, nrows
      character(len=STRING_SIZE), dimension(:), allocatable :: output

      character(len=LINE_SIZE) :: line
      integer :: alloc_stat, io_stat, i

      allocate(output(nrows), stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation error in <in <read_list>"

      call skip_line(eunit)
      do i = 1, nrows
         line = get_line(eunit)

         read(line, *, iostat=io_stat) output(i)
         if (io_stat /= 0) then
            print *, '>>> ' // adjustl(trim(line))
            stop 'read error in <get_list>'
         end if
      end do
      call skip_line(eunit)
   end function read_list

   !-----------------------------------
   ! FUNCTION: read_matrix
   !-----------------------------------
   function read_matrix(eunit, nrows, ncols) result(output)
      integer, intent(in) :: eunit
      integer, intent(in) :: nrows
      integer, intent(in) :: ncols
      real(8), dimension(:, :), allocatable :: output

      character(len=LINE_SIZE) :: line
      integer :: alloc_stat, io_stat
      integer :: i, j

      allocate(output(nrows, ncols), stat=alloc_stat)
      if (alloc_stat /= 0) stop "allocation error in <in <read_matrix>"

      call skip_line(eunit)
      do i = 1, nrows
         call skip_line(eunit)
         do j = 1, ncols
            line = get_line(eunit)

            read(line, *, iostat=io_stat) output(i, j)
            if (io_stat /= 0) then
               print *, '>>> ' // adjustl(trim(line))
               stop 'read error in <read_matrix>'
            end if
         end do
         call skip_line(eunit)
      end do
      call skip_line(eunit)
   end function read_matrix

   !-----------------------------------
   ! FUNCTION: read_real
   !-----------------------------------
   function read_real(eunit) result(output)
      integer, intent(in) :: eunit
      real(8) :: output

      character(len=LINE_SIZE) :: line
      integer :: io_stat
      integer :: left

      line = get_line(eunit)
      left = scan(line, ':') + 1

      read(line(left:), *, iostat=io_stat) output
      if (io_stat /= 0) then
         print *, '>>> ' // adjustl(trim(line))
         stop 'read error in <read_real>'
      end if
   end function read_real

   !-----------------------------------
   ! FUNCTION: read_string
   !-----------------------------------
   function read_string(eunit) result(output)
      integer, intent(in) :: eunit
      character(len=STRING_SIZE) :: output

      character(len=LINE_SIZE) :: line
      integer :: io_stat, left

      line  = get_line(eunit)
      left  = scan(line, ':') + 1

      read(line(left:), *, iostat=io_stat) output
      if (io_stat /= 0) then
         print *, '>>> ' // adjustl(trim(line))
         stop 'read error in <get_string>'
      end if
   end function read_string

   !-----------------------------------
   ! FUNCTION: read_vector
   !-----------------------------------
   function read_vector(eunit, nrows) result(output)
      integer, intent(in) :: eunit, nrows
      real(8), dimension(:), allocatable :: output

      character(len=LINE_SIZE) :: line
      integer :: alloc_stat, io_stat
      integer :: i

      allocate(output(nrows), stat=alloc_stat)
      if (alloc_stat /= 0) stop 'allocation error in <read_vector>'

      call skip_line(eunit)
      do i = 1, nrows
         line = get_line(eunit)

         read(line, *, iostat=io_stat) output(i)
         if (io_stat /= 0) then
            print *, '>>> ' // adjustl(trim(line))
            stop 'read error in <read_vector>'
         end if
      end do
      call skip_line(eunit)
   end function read_vector

   !-----------------------------------
   ! SUBROUTINE: skip_line
   !-----------------------------------
   subroutine skip_line(eunit)
      integer, intent(in) :: eunit

      character(len=LINE_SIZE) :: line

      line = get_line(eunit)
   end subroutine skip_line

end module read_json_module
