!--------------------------------------------------------------------

!--------------------------------------------------------------------

module basekinds

!--------------------------------------------------------------------

!--------------------------------------------------------------------

!

!--------------------------------------------------------------------

!< Description:

!< This module defines the parameters used 
!< intended to match most blas/LApack definition
!< of single and double precision.
!< and a reasonably sized kind parameter for integers.

!< Also, machine precision is defined twice, for double
!< and single precision

!< Lastly, there are also subroutines to define
!< file unit numbers and to quicksort 1D array indexes

!--------------------------------------------------------------------

!

!--------------------------------------------------------------------

! Modules and Global Variables

!--------------------------------------------------------------------

!--------------------------------------------------------------------

!

  implicit none

!

!--------------------------------------------------------------------

!--------------------------------------------------------------------

! Fixed Kind parameters

!--------------------------------------------------------------------

!! single precision parameter

  integer, parameter :: kind_single = kind(16e0)

!! double precision parameter

  integer, parameter :: &

  & kind_double = selected_real_kind(2*precision(16.0_kind_single))

!! 4 byte parameter for integers (int)

  integer, parameter :: kind_integer = 4

!! real, double precision epsilon (machine precision)

  real(kind_double), parameter :: &

  & eps_double = epsilon(real(0,kind=kind_double))

!! real, single precision epsilon (machine precision)

  real(kind_double), parameter :: &

  & eps_single = epsilon(real(0,kind=kind_single))


!! real, double precision log10(epsilon) (machine precision)

  real(kind_double), parameter :: &
  & logeps_double = log10(epsilon(real(0,kind=kind_double)))

!! real, double precision pi

  real(kind_double), parameter :: &
  & pi = atan(real(1,kind=kind_double))*4

!! bohr converted to double precision, note: accuracy at best single precision
  real(kind_double), parameter :: &
  & bohr = 5.29177210903d-1

!--------------------------------------------------------------------

contains

!-------------------------------------------------------------------- 

  subroutine find_free_file_unit(funit,ierr) 

!-------------------------------------------------------------------- 

! 

!-------------------------------------------------------------------- 

!< Description: 

!< This routine returns a free file unit, 14 < funit < 256 

!-------------------------------------------------------------------- 

! 

    implicit none 

! 

!-------------------------------------------------------------------- 

! Output Parameters 

!-------------------------------------------------------------------- 

!! name of file to be read 

    integer(kind_integer), intent(out) :: funit 

!-------------------------------------------------------------------- 

! Error Parameter 

!-------------------------------------------------------------------- 

    integer(kind_integer), intent(inout) :: ierr 

!-------------------------------------------------------------------- 

!  Local Variables 

!-------------------------------------------------------------------- 

!! dummy indexes 

    integer(kind_integer) :: j = 0 

!! logic for file unit 

    logical :: already_used = .true. 

!-------------------------------------------------------------------- 

 

! Assume file unit number 0 to 14 are reserved for other output files

!! INCREASE THIS NUMBER TO RESERVE MORE FILE UNIT NUMBERS

    j = 14 

    already_used = .true. 

    ierr = -1

    do  

      j = j + 1 

      inquire(unit=j,opened=already_used) 

      !! Example of how to skip a file unit number if you so choose

     ! if (j.eq.10) cycle 

      !!! automatically stop if more than 256 funits are used. 

      !!! Think carefully before increasing this number

      if (j.gt.256) exit

      if(.not.already_used) then 

        funit = j

        ierr = 0 

        exit 

      end if 

    end do 

 

!-------------------------------------------------------------------- 

  end subroutine find_free_file_unit 

!-------------------------------------------------------------------- 

!--------------------------------------------------------------------

!--------------------------------------------------------------------

  recursive subroutine quicksort_rdp(n,obj,dex,first,last)

!--------------------------------------------------------------------

!--------------------------------------------------------------------

!

!--------------------------------------------------------------------

! Description:

!--------------------------------------------------------------------

!! reference to: https://gist.github.com/t-nissie/479f0f16966925fa29ea

!! sorts a linear array of real(kind_double), obj, 

!! from smallest to largest

!! dex contains the ordering, to match the new order to the old,

!! i.e. dex(1) is the original position of the smallest value

!--------------------------------------------------------------------

!

!--------------------------------------------------------------------
! Modules and Global Variables
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!
    implicit none
!
!--------------------------------------------------------------------

! Input/Output Parameters

!--------------------------------------------------------------------

!! size of obj and dex

    integer(kind_integer), intent(in) :: n

!!  array to be sorted

    real(kind_double), intent(inout) :: obj(n)

!! indexing of array to be sorted

    integer(kind_integer), intent(inout) :: dex(n)

!--------------------------------------------------------------------

! Input Parameters

!--------------------------------------------------------------------

!! first element

    integer(kind_integer), intent(in) :: first

!! the last element

    integer(kind_integer), intent(in) :: last
!--------------------------------------------------------------------

!  Local Variables

!--------------------------------------------------------------------

!!  pivot value, by default is the value in the middle of first and last

    real(kind_double) :: p

!!  dummy variable for copying

    real(kind_double) :: t

!!  dummy variable for copying

    integer(kind_integer) :: k = 0

!!  integer for do loops

    integer(kind_integer) :: j1,j2 = 0

!--------------------------------------------------------------------

    p = obj(int((first+last)*0.5,kind=kind_integer))

    j1 = first

    j2 = last

    do

      do while (obj(j1).lt.p)

        j1 = j1+1

      end do

      do while (obj(j2).gt.p)

        j2 = j2-1

      end do

      if (j1.ge.j2) exit

      t = obj(j1)

      k = dex(j1)

      obj(j1) = obj(j2)

      dex(j1) = dex(j2)

      obj(j2) = t

      dex(j2) = k

      j1 = j1+1

      j2 = j2-1

    end do

    if (first.lt.j1-1) then

      call quicksort_rdp(n,obj,dex,first,j1-1)

    end if

    if (last.gt.j2+1) then

      call quicksort_rdp(n,obj,dex,j2+1,last)

    end if

!--------------------------------------------------------------------
!--------------------------------------------------------------------

  end subroutine quicksort_rdp

!--------------------------------------------------------------------
!--------------------------------------------------------------------




!--------------------------------------------------------------------
!--------------------------------------------------------------------

end module basekinds

!--------------------------------------------------------------------
!--------------------------------------------------------------------
