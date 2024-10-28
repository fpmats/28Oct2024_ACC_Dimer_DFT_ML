!--------------------------------------------------------------------
!--------------------------------------------------------------------
program near_gap_number
!--------------------------------------------------------------------
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! Modules and Global Variables
!--------------------------------------------------------------------
use basekinds
!--------------------------------------------------------------------
!
  implicit none
!
!--------------------------------------------------------------------

!! error variable
  integer(kind_integer) ::ierr = 0

!! loop variables
  integer(kind_integer) :: j,k,l,m = 0
  integer(kind_integer) :: a,b,c,d = 0

!! user options
  character(len=32) :: input = ''
  character(len=32) :: input2 = ''
  character(len=32) :: file_name = 'orbital_energies.out'
  real(kind_double) :: window = 7.34987e-3
  integer(kind_integer) :: total_mo = 0
  integer(kind_integer) :: homo = 0
  integer(kind_integer) :: lines = 0
  integer(kind_integer) :: remainder = 0
  logical :: occupied = .false.

!! file opening
  logical :: file_exists = .false.
  integer(kind_integer) :: funit = 0
  character(len=80) :: file_string = ''
  character(len=16) :: word_string(10) = ''
  integer(kind_integer) :: n_words = 0
  real(kind_double), allocatable :: eigenvalue(:)

!! checking command line options:
  j = command_argument_count()
!! loop over command line
  k = 1

  if (j.gt.0) then

    do 
      call get_command_argument(k,value=input,status=ierr)
      if (ierr.ne.0) stop

      if ((input.eq.'-help').or.(input.eq.'--help')) then

        print *, 'Program which returns the '
        print *, 'number of near-gap KS orbitals based on their eigenvalues'
        print *, 'from a g16.fchk eigenvalue data block on file'
        print *, ''
        print *, 'Maximum file name length is 80 (adjusted in source code)'
        print *, ''

        print *, 'options:'
        print *, '--help             display this message'
        print *, ''

        print *, '-file <*.out>      file name for file containing'
        print *, '                   the block of data.'
        print *, '                   Default: <*.out> = orbital_energies.out'
        print *, '                   Limit of 32 character filename.'
        print *, ''

        print *, '-homo <m>          index of HOMO. *REQUIRED*'
        print *, '                   <m> must be a positive integer'
        print *, '                   less than number of eigenvalues.'
        print *, ''

        print *, '-occupied          indicates that number of occupied'
        print *, '                   orbitals in the window from the'
        print *, '                   HOMO and down is returned.'
        print *, '                   Without this line, number of'
        print *, '                   unoccupied orbitals in the window'
        print *, '                   from the LUMO and up is returned.'
        print *, '                   Number returned does not include this'
        print *, '                   reference orbital.'
        print *, ''

        print *, '-window <x>        changes the size of the window'
        print *, '                   used. <x> must be a real number.'
        print *, '                   units are in a.u. , Hartrees.'
        print *, '                   Default: <x> = 7.34987e-3'
        print *, '                   corresponds to about 0.2eV'
        print *, ''

        stop

      else if (input.eq.'-file') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) file_name
        if (ierr.ne.0) exit

      else if (input.eq.'-homo') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) homo
        if (ierr.ne.0) exit

      else if (input.eq.'-occupied') then

        occupied = .true.

      else if (input.eq.'-window') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) window
        if (ierr.ne.0) exit

      else if (input.eq.'>') then
        exit
      else if (input.eq.'>>') then
        exit
      else if (input.eq.'|') then
        exit
      end if

      k = k + 1
      if (k.gt.j) exit

    end do
  end if 

!! sanity check
  if (ierr.ne.0) then
    print *, 'faulty command line input failed to match options!'
    stop
  end if
  if (file_name.eq.'') then
    print *, 'No file name provided!'
    stop
  end if
  if (homo.eq.0) then
    print *, 'no HOMO specified!'
    stop
  end if


!! open input file
  call find_free_file_unit(funit,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units!'
    stop
  end if
  inquire(file=file_name,exist=file_exists)
  if (.not.file_exists) then
    print *, 'file to be read does not exist!'
    stop
  end if
  open(unit=funit,file=file_name,action='read',iostat=ierr)

!! read info on data block
  read(unit=funit,fmt='(a)',iostat=ierr) file_string
  if (ierr.ne.0) then
    print *, 'Problem reading info line'
    close(unit=funit,iostat=ierr,status='keep')
    stop
  end if
  call parse_line(file_string,word_string,n_words)
  if (n_words.ne.6) then
    print *, 'Weird number of words in info line'
    close(unit=funit,iostat=ierr,status='keep')
    stop
  end if
  read(word_string(6),*,iostat=ierr) total_mo

!! allocate eigenvalues
  allocate(eigenvalue(total_mo))

!! number of lines of data block, rounded down if needed
  lines = total_mo/5
!! number of elements in last incomplete line, if any
  remainder = modulo(total_mo,5)

!! read in values to eigenvalues
  l = 1
  do j = 1, lines
    read(unit=funit,fmt='(a)',iostat=ierr) file_string
    if (ierr.ne.0) then
      print *, 'Problem reading data line'
      close(unit=funit,iostat=ierr,status='keep')
      deallocate(eigenvalue)
      stop
    end if
    call parse_line(file_string,word_string,n_words)
    if (n_words.ne.5) then
      print *, 'Weird number of words in line ',j+1
      close(unit=funit,iostat=ierr,status='keep')
      deallocate(eigenvalue)
      stop
    end if
    do k = 1, 5
      read(word_string(k),*,iostat=ierr) eigenvalue(l) 
      l = l + 1
    end do
  end do
! read incomplete last line, if it exists
  if (remainder.ne.0) then
    read(unit=funit,fmt='(a)',iostat=ierr) file_string
    if (ierr.ne.0) then
      print *, 'Problem reading last line'
      close(unit=funit,iostat=ierr,status='keep')
      deallocate(eigenvalue)
      stop
    end if
    call parse_line(file_string,word_string,n_words)
    if (n_words.ne.remainder) then
      print *, 'Weird number of words in last line'
      close(unit=funit,iostat=ierr,status='keep')
      deallocate(eigenvalue)
      stop
    end if
    do k = 1, remainder
      read(word_string(k),*,iostat=ierr) eigenvalue(l) 
      l = l + 1
    end do
  end if

!! check all values are read
  if (total_mo.ne.(l-1)) then
    print *, 'Number of values in file not consistent'
    close(unit=funit,iostat=ierr,status='keep')
    deallocate(eigenvalue)
    stop
  end if

  if (window.gt.0) then
    if (occupied) then
  !! occupied analysis
      do j = homo-1, 1, -1
        if ((eigenvalue(homo) - eigenvalue(j) - window).gt.1.0d-9) then
          print *, homo-j-1
          exit
        end if
      end do
    else
  !! unoccupied analysis
      do j = homo+2, total_mo
        if ((eigenvalue(j) - eigenvalue(homo+1) - window).gt.1.0d-9) then
          print *, j-homo-2
          exit
        end if
      end do
    end if
  else
    if (occupied) then
  !! occupied analysis
      do j = homo+1, total_mo
        if ((eigenvalue(j) - eigenvalue(homo) + window).gt.1.0d-9) then
          print *, j-homo-1
          exit
        end if
      end do
    else
  !! unoccupied analysis
      do j = homo, 1, -1
        if ((eigenvalue(homo+1) - eigenvalue(j) + window).gt.1.0d-9) then
          print *, homo-j
          exit
        end if
      end do
    end if
  end if 

!! close file and deallocate at end of program
  close(unit=funit,iostat=ierr,status='keep')
  deallocate(eigenvalue)

!--------------------------------------------------------------------
!--------------------------------------------------------------------
end program near_gap_number
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!! subroutine to separate words in file string into an array
!! dimensions are fixed but easily changed
!! single character separator defined internally
subroutine parse_line(file_string,word_string,n_words)

use basekinds
  implicit none

  character(len=80), intent(in) :: file_string
  character(len=16), intent(out) :: word_string(10)
  integer(kind_integer), intent(out) :: n_words

  integer(kind_integer) :: j, k, l
  character(len=1) :: separator = ' '

  l = 1
  j = 0
  do
    if (j.eq.len(file_string)) exit
    j = j + 1
    if (file_string(j:j).ne.separator) then
      k = j
      if (j.eq.len(file_string)) then
        word_string(l) = file_string(j:j)
        l = l + 1
        exit
      end if
      do
        j = j + 1
        if (file_string(j:j).eq.separator) then
          word_string(l) = file_string(k:j-1)
          l = l + 1
          exit
        else if (j.eq.len(file_string)) then
          word_string(l) = file_string(k:j)
          l = l + 1
          exit
          exit
        end if
      end do
    end if
  end do 

  n_words = l - 1

end subroutine parse_line

subroutine element_to_lower(strIn,strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page
use basekinds

  implicit none

  character(len=2), intent(in) :: strIn
  character(len=2) :: strOut
  integer(kind_integer) :: i,j

  do i = 1, 2
    j = iachar(strIn(i:i))
    if (j.ge.iachar("A") .and. j.le.iachar("Z") ) then
      strOut(i:i) = achar(iachar(strIn(i:i))+32)
    else
      strOut(i:i) = strIn(i:i)
    end if
  end do

end subroutine element_to_lower


