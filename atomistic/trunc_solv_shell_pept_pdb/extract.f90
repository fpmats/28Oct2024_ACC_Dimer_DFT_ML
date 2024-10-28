!--------------------------------------------------------------------
!--------------------------------------------------------------------
program extract
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
  integer(kind_integer) :: j,k,l,m,n = 0
  integer(kind_integer) :: a,b = 0

!! user options
  character(len=32) :: input = ''
  character(len=32) :: input2 = ''
  character(len=32) :: file_name = ''
  character(len=64) :: temp_name = ''
  character(len=64) :: temp2_name = ''
  character(len=64) :: out_name = ''
  character(len=64) :: coord_name = ''
  character(len=64) :: xyz_name = ''
  character(len=64) :: com_name = ''
  integer(kind_integer) :: header_lines = 5

!! file opening
  integer(kind_integer) :: funit,funit2,funit3 = 0
  integer(kind_integer) :: funit4,funit5,funit6 = 0
  integer(kind_integer) :: funit7 = 0
  logical :: file_exists = .false.
  character(len=80) :: file_string = ''
  character(len=80) :: file_string2 = ''
  character(len=80) :: file_string3 = ''
  character(len=80) :: file_string4 = ''
  character(len=8) :: word_string(13) = ''
  integer(kind_integer) :: n_words = 0

  integer(kind_integer) :: n_atoms = 0
  integer(kind_integer) :: n_peptide_atoms = 0
  integer(kind_integer) :: n_shell_atoms = 0
  character(len=8) :: str_peptide_atoms = ''
  integer(kind_integer) :: n_charge = 0
  character(len=8) :: str_charge = ''
  character(len=2) :: str_element = ''

!! atom of interest
  logical :: atom_of_interest = .false.
  logical :: selective = .false.
  logical :: selective_charge = .false.
  logical :: selective_hetero = .false.

  real(kind_double), allocatable :: peptide_pos(:,:) 
  character(len=2) , allocatable :: peptide_element(:) 
  real(kind_double) :: atom_pos(3) = 0.0d0
  real(kind_double) :: atom_pos_bohr(3) = 0.0d0
  real(kind_double) :: distance_vec(3) = 0.0d0
  real(kind_double) :: distance = 0.0d0
  real(kind_double) :: cutoff = 2.90d0
  real(kind_double) :: threshold = 0.0d0
  real(kind_double) :: cutoff_N = 3.14d0
  real(kind_double) :: cutoff_O = 3.35d0
  real(kind_double) :: cutoff_S = 4.25d0


!! checking command line options:
  j = command_argument_count()
!! loop over command line
  k = 1

  if (j.gt.0) then

    do 
      call get_command_argument(k,value=input,status=ierr)
      if (ierr.ne.0) stop

      if ((input.eq.'-help').or.(input.eq.'--help')) then

        print *, 'Truncating Solvation Shell of a peptide pdb file,' 
        print *, ' from MD in aqueous solvent with counter ions.'
        print *, '  outputs, reads and deletes a pdb file with just'
        print *, '   the peptide atoms, also outputs a pdb file '
        print *, '    with all desired components:'
        print *, '     the peptide, first solvent shell with '
        print *, '      both water and counter ions.'
        print *, ''
        print *, 'An .xyz file with only the peptide is produced.'
        print *, ''
        print *, 'A Gaussian16 input file with the peptide and only'
        print *, ' the water solvent shell is produced as well.'
        print *, ''

        print *, 'options (in the order not to overwrite each other):'
        print *, '--help             display this message'
        print *, ''

        print *, '-file <file>       file name for input file containing'
        print *, '                   the initial pdb data. *REQUIRED*'
        print *, '                   Limit of 32 character string.'
        print *, ''

        print *, '-temp <temp>       file name for temp file to'
        print *, '                   contain only the peptide pdb data.'
        print *, '                   Default: "temp_"//"<file>"'
        print *, '                   Limit of 64 character string.'
        print *, '                   Will be overwritten if it exists.'
        print *, ''

        print *, '-temp2 <temp2>     file name for temp file to'
        print *, '                   contain only the interesting pdb data.'
        print *, '                   Default: "temp2_"//"<file>"'
        print *, '                   Limit of 64 character string.'
        print *, '                   Will be overwritten if it exists.'
        print *, ''

        print *, '-out <out>         file name for output file to'
        print *, '                   contain desired output pdb data.'
        print *, '                   Default: "out_"//"<file>"'
        print *, '                   Limit of 64 character string.'
        print *, '                   Will be overwritten if it exists.'
        print *, ''

        print *, '-coord <coord>     file name for output coord file to'
        print *, '                   contain desired output data in bohrs.'
        print *, '                   Default: "<file>"//".coord"'
        print *, '                   Limit of 64 character string.'
        print *, '                   Will be overwritten if it exists.'
        print *, ''

        print *, '-xyz <xyz>         file name for output xyz file to'
        print *, '                   contain peptide only xyz data.'
        print *, '                   Default: "<file>"//".xyz"'
        print *, '                   Limit of 64 character string.'
        print *, '                   Will be overwritten if it exists.'
        print *, ''

        print *, '-com <com>         file name for output com file to'
        print *, '                   contain Gaussian16 input file.'
        print *, '                   Default: "<file>"//".com"'
        print *, '                   Limit of 64 character string.'
        print *, '                   Will be overwritten if it exists.'
        print *, ''

        print *, '-cutoff <cut>      real number for threshold for'
        print *, '                   environment molecule to be'
        print *, '                   considered part of the first'
        print *, '                   solvation shell of peptide.'
        print *, '                   Based on oxygen atom of water,'
        print *, '                   and central atom of counter-ion.'
        print *, '                   (Complex ions must individually'
        print *, '                    added to the source code.)'
        print *, '                   This overrides any special cutoff'
        print *, '                   for special atoms set before.'
        print *, '                   Default: 2.9d0 Angstroms'
        print *, ''

        print *, '-header <head>     integer for number of header lines'
        print *, '                   in the input file,'
        print *, '                   which must be skipped to start'
        print *, '                   reading coordinates.'
        print *, '                   Default: 5'
        print *, ''

        print *, '-selective_charge  by default a uniform solvation '
        print *, '                   shell around all peptide atoms is created,'
        print *, '                   including this option only considers'
        print *, '                   hetero-atoms O,N of charged functional'
        print *, '                   groups LYS and GLU, and peptide termini.'
        print *, '                   See following options for different cutoffs'
        print *, '                   for different peptide atoms.'
        print *, '                   Mutually exclusive from other '
        print *, '                   selective options as the groups are'
        print *, '                   not exclusive.'
        print *, ''

        print *, '-selective_hetero  by default a uniform solvation '
        print *, '                   shell around all peptide atoms is created,'
        print *, '                   including this option only considers'
        print *, '                   hetero-atoms O,N,S.'
        print *, '                   See following options for different cutoffs'
        print *, '                   for different peptide atoms.'
        print *, '                   Mutually exclusive from other '
        print *, '                   selective options as the groups are'
        print *, '                   not exclusive.'
        print *, ''

        print *, '-cut_N <cut_N>     real number for threshold for'
        print *, '                   environment molecule to be'
        print *, '                   considered part of the first'
        print *, '                   solvation shell of selected N'
        print *, '                   (described in -selective option).'
        print *, '                   This is overwritten by any -cutoff'
        print *, '                   set after this option.'
        print *, '                   Default: 3.14d0 Angstroms'
        print *, ''

        print *, '-cut_O <cut_O>     real number for threshold for'
        print *, '                   environment molecule to be'
        print *, '                   considered part of the first'
        print *, '                   solvation shell of selected O'
        print *, '                   (described in -selective option).'
        print *, '                   This is overwritten by any -cutoff'
        print *, '                   set after this option.'
        print *, '                   Default: 3.35d0 Angstroms'
        print *, ''

        print *, '-cut_S <cut_S>     real number for threshold for'
        print *, '                   environment molecule to be'
        print *, '                   considered part of the first'
        print *, '                   solvation shell of selected S'
        print *, '                   (described in -selective option).'
        print *, '                   This is overwritten by any -cutoff'
        print *, '                   set after this option.'
        print *, '                   Default: 4.25d0 Angstroms'
        print *, ''


        stop

      else if (input.eq.'-file') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) file_name
        print *, 'Reading pdb data in file: ',input2
        if (ierr.ne.0) exit
        temp_name = 'temp_'//trim(file_name)
        temp2_name = 'temp2_'//trim(file_name)
        out_name = 'out_'//trim(file_name)
        coord_name = trim(file_name)//'.coord'
        xyz_name = trim(file_name)//'.xyz'
        com_name = trim(file_name)//'.com'

      else if (input.eq.'-temp') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) temp_name
        print *, 'Temp peptide pdb data written to file: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-temp2') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) temp2_name
        print *, 'Temp interesting pdb data written to file: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-out') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) out_name
        print *, 'Output pdb data written to file: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-coord') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) coord_name
        print *, 'Output coord data written to file: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-xyz') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) xyz_name
        print *, 'Output peptide xyz data written to file: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-com') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) com_name
        print *, 'Gaussian16 input written to file: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-cutoff') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) cutoff
        cutoff_N = cutoff
        cutoff_O = cutoff
        cutoff_S = cutoff
        print *, 'Cutoff for first solvation shell is: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-cut_N') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) cutoff_N
        print *, 'Nitrogen cutoff for first solvation shell is: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-cut_O') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) cutoff_O
        print *, 'Oxygen cutoff for first solvation shell is: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-cut_S') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) cutoff_S
        print *, 'Sulfur cutoff for first solvation shell is: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-header') then

        k = k + 1
        call get_command_argument(k,value=input2,status=ierr)
        if (ierr.ne.0) exit
        read(input2,*,iostat=ierr) header_lines
        print *, 'Number of header lines in input file is: ',input2
        if (ierr.ne.0) exit

      else if (input.eq.'-selective_charge') then

        selective_charge = .true.
        print *, 'Selective solvation shell for charged groups constructed.'
        if (ierr.ne.0) exit

      else if (input.eq.'-selective_hetero') then

        selective_hetero = .true.
        print *, 'Selective solvation shell for hetero-atoms constructed.'
        if (ierr.ne.0) exit

      else if (input.eq.'>') then
        exit
      else if (input.eq.'>>') then
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

  if (selective_charge.and.selective_hetero) then
    print *, 'Two exclusive -selective options selected!'
    stop
  else if (selective_charge.or.selective_hetero) then
    selective = .true.
  end if


!! open input file
  call find_free_file_unit(funit,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units! (1)'
    stop
  end if
  inquire(file=file_name,exist=file_exists)
  if (.not.file_exists) then
    print *, 'file to be read does not exist!'
    stop
  end if
  open(unit=funit,file=file_name,action='read',iostat=ierr)



!! open temp file, copy unit cell info from input
  call find_free_file_unit(funit2,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units! (2)'
    stop
  end if
  open(unit=funit2,file=temp_name,action='write',iostat=ierr)
!! open temp2 file, copy unit cell info from input
  call find_free_file_unit(funit6,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units! (6)'
    stop
  end if
  open(unit=funit6,file=temp2_name,action='write',iostat=ierr)
!! open output file, copy unit cell info from input
  call find_free_file_unit(funit3,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units! (3)'
    stop
  end if
  open(unit=funit3,file=out_name,action='write',iostat=ierr)
!! open output coord file.
  call find_free_file_unit(funit7,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units! (7)'
    stop
  end if
  open(unit=funit7,file=coord_name,action='write',iostat=ierr)
!! open output xyz file.
  call find_free_file_unit(funit4,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units! (4)'
    stop
  end if
  open(unit=funit4,file=xyz_name,action='write',iostat=ierr)
!! open output com file.
  call find_free_file_unit(funit5,ierr)
  if (ierr.ne.0) then
    print *, 'No free file units! (5)'
    stop
  end if
  open(unit=funit5,file=com_name,action='write',iostat=ierr)

!! read and repeat cell info
  do j = 1, header_lines
    read(unit=funit,fmt='(a)',iostat=ierr) file_string
    if (ierr.ne.0) then
      print *, 'Problem reading cell info, line:', j
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
    if (j.eq.3) then
      write(unit=funit2,fmt='(a)',iostat=ierr) &
  &     'REMARK    THIS IS A PEPTIDE SUBSYSTEM FROM A SIMULATION BOX'
      write(unit=funit6,fmt='(a)',iostat=ierr) &
  &     'REMARK    THIS IS INTERESTING ATOMS OF PEPTIDE FROM A SIMULATION BOX'
      write(unit=funit3,fmt='(a)',iostat=ierr) &
  &     'REMARK    THIS IS A PEPTIDE AND SHELL SUBSYSTEM FROM A SIMULATION BOX'
    else
      write(unit=funit2,fmt='(a)',iostat=ierr) file_string
      write(unit=funit6,fmt='(a)',iostat=ierr) file_string
      write(unit=funit3,fmt='(a)',iostat=ierr) file_string
    end if
  end do

!! read and write peptide list to temp
  j = 0
  k = 0
  n = 0
  a = 0
  b = 0 !! number of atoms associated with this charged residue
!! n_charge = 0 based on initialization, presume
!! charge on both ends of peptide cancel.
!! One needs to manually include checks for modified peptide ends
!! which change the 'default' zwitter-ionic nature
!! of the peptide backbone
  do
    j = j + 1
    read(unit=funit,fmt='(a)',iostat=ierr) file_string
    if (ierr.ne.0) then
      print *, 'Problem reading atom index: ',j
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
    call parse_PDB_line(file_string,word_string,n_words)
    if (n_words.le.0) then
      print *, 'Problem reading atom index: ',j
      print *, file_string
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
    if (trim(word_string(1)).eq.'END') then
      print *, 'End of file at line: ',j
      print *, 'No solvent in file.'
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
    if (trim(word_string(1)).eq.'TER') then !! TODO, may not be termination line!
      print *, 'Early termination line of file: ',j
    end if
!! Assuming all peptide chains have a chain ID
    if (trim(word_string(5)).eq.'') then ! Not a peptide atom 
      if (selective_charge .and. (a.eq.-4)) then !! passed over terminal COO group and off the peptide
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string2
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string3
        n = n + 2
      end if
      exit
    end if
!! if not in a charged residue sequence, check if residue is charged
!! add charge to total, start counting atoms to decide when to check
!! condition again.
! catch for COO terminal charged residue
    if (a.eq.-5) then
      if (trim(word_string(3)).eq.'OT2') then !! terminal? O
        n_charge = n_charge - 1
      else
        a = 0 !! restart checks for charged residue in this loop
      end if
    else if (a.eq.-4) then  !! passed over terminal COO(?) group
      if (trim(word_string(3)).eq.'HT1') then !! terminal H
        n_charge = n_charge + 1 !! cancel out pre-emptive charge addition
      else
        if (selective_charge) then
          write(unit=funit6,fmt='(a)',iostat=ierr) file_string2
          write(unit=funit6,fmt='(a)',iostat=ierr) file_string3
          n = n + 2
        end if
        a = 0 !! restart checks for charged residue in this loop
      end if
    else if (a.eq.-3) then  !! passed over terminal COOH group
      a = 0 !! restart checks for charged residue in this loop
    end if
    if (a.eq.0) then
      if (trim(word_string(4)).eq.'ASP') then
        n_charge = n_charge - 1
        a = 1
        b = 12
      else if (trim(word_string(4)).eq.'GLU') then
        n_charge = n_charge - 1
        a = 1
        b = 15
      else if (trim(word_string(4)).eq.'HIS') then
        n_charge = n_charge + 1
        a = 1
        b = 15
      else if (trim(word_string(4)).eq.'LYS') then
        n_charge = n_charge + 1
        a = 1
        b = 22
      else if (trim(word_string(4)).eq.'ARG') then
        n_charge = n_charge + 1
        a = 1
        b = 24
      else if (trim(word_string(4)).eq.'LEU') then
      else if (trim(word_string(4)).eq.'ALA') then
      else if (trim(word_string(4)).eq.'ILE') then
      else if (trim(word_string(4)).eq.'GLN') then
      else if (trim(word_string(4)).eq.'PHE') then
      else
        print *, 'Unrecognized Residue: ', trim(word_string(4))
        stop
      end if
    else
      if (trim(word_string(3)).eq.'HT2') then
        b = b + 1 !! add (H) for N-terminus residue
      else if (trim(word_string(3)).eq.'HT3') then
        if (selective_charge) then
          write(unit=funit6,fmt='(a)',iostat=ierr) file_string4
          n = n + 1
        end if
        b = b + 1 !! add (H) for N-terminus residue
        n_charge = n_charge + 1 !! and associated charge
      end if
      a = a + 1
      if (a.eq.b) then
        a = -5
        b = 0
      end if
    end if
! at this point the file line is confirmed to be part of the peptide,
! write line to temp and output files.
    write(unit=funit2,fmt='(a)',iostat=ierr) file_string
    write(unit=funit3,fmt='(a)',iostat=ierr) file_string
!! LOGIC FOR SELECTING ATOMS OF INTEREST
!! ALL CHARGED
    if (selective_charge) then
      if ((trim(word_string(3)).eq.'NZ').and.(trim(word_string(4)).eq.'LYS')) then
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
        n = n + 1
      else if ((trim(word_string(3)).eq.'OE1').and.(trim(word_string(4)).eq.'GLU')) then
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
        n = n + 1
      else if ((trim(word_string(3)).eq.'OE2').and.(trim(word_string(4)).eq.'GLU')) then
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
        n = n + 1
!!      else if ((trim(word_string(3)).eq.'OE1').and.(trim(word_string(6)).eq.'28')) then
!!        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
!!        n = n + 1
!!      else if ((trim(word_string(3)).eq.'NE2').and.(trim(word_string(6)).eq.'28')) then
!!        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
!!        n = n + 1
      end if
    else if (selective_hetero) then
!! ALL HETERO
      if (trim(word_string(13)).eq.'N') then
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
        n = n + 1
      else if (trim(word_string(13)).eq.'O') then
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
        n = n + 1
      else if (trim(word_string(13)).eq.'S') then
        write(unit=funit6,fmt='(a)',iostat=ierr) file_string
        n = n + 1
      end if
    end if
    file_string4 = file_string3
    file_string3 = file_string2
    file_string2 = file_string
  end do

!! peptide temp file is now complete, can be closed and reopened
!! as read only.
  write(unit=funit2,fmt='(a)',iostat=ierr) 'END'
  close(unit=funit2,iostat=ierr,status='keep')
!! TODO: Improve funit 6 opening and closing
  write(unit=funit6,fmt='(a)',iostat=ierr) 'END'
  close(unit=funit6,iostat=ierr,status='keep')
!! The number of peptide atoms is now determined.
  n_peptide_atoms = j - 1
  print *, 'Number of peptide atoms: ', n_peptide_atoms
  print *, 'Charge of peptide: ', n_charge

!! We can now allocate arrays for the peptide positions
!!  to be used for scanning.
  print *, 'Atoms to be solvated by `-selective_*` option: ', n
  if (.not.selective) n = n_peptide_atoms
  print *, 'Atoms to be solvated: ', n
  allocate(peptide_pos(3,n)) 
  allocate(peptide_element(n)) 

!! Block to read in temp peptide data and print xyz and com file.
  open(unit=funit2,file=temp_name,action='read',iostat=ierr)

  write(str_peptide_atoms,'(i8)') n_peptide_atoms
  write(str_charge,'(i8)') n_charge
  write(unit=funit4,fmt='(a)',iostat=ierr) adjustl(str_peptide_atoms)
  write(unit=funit4,fmt='(a)',iostat=ierr) ''

  write(unit=funit5,fmt='(a)',iostat=ierr) '#P PBE1PBE/cc-pVDZ int(ultrafine) SCRF(CPCM) NoSymm'
  write(unit=funit5,fmt='(a)',iostat=ierr) 'Pop=(Orbitals=20,ThreshOrbitals=1)'
  write(unit=funit5,fmt='(a)',iostat=ierr) ''
  write(unit=funit5,fmt='(a)',iostat=ierr) 'generated from: '//file_name 
  write(unit=funit5,fmt='(a)',iostat=ierr) ''
  write(unit=funit5,fmt='(a)',iostat=ierr) trim(adjustl(str_charge))//' 1'

  write(unit=funit7,fmt='(a)',iostat=ierr) '$coord'
  
  do l = 1, header_lines
    read(unit=funit2,fmt='(a)',iostat=ierr) file_string2
    if (ierr.ne.0) then
      print *, 'Problem reading temp file info, line: ', l
      deallocate(peptide_pos) 
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
  end do
  do k = 1, n_peptide_atoms
    read(unit=funit2,fmt='(a)',iostat=ierr) file_string2
    if (ierr.ne.0) then
      print *, 'Problem reading temp peptide atom index: ', k
      deallocate(peptide_pos) 
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
    call parse_PDB_line(file_string2,word_string,n_words)
    if (n_words.le.0) then
      print *, 'Problem reading temp peptide atom index: ', k
      print *, file_string2
      deallocate(peptide_pos) 
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
    if (trim(word_string(1)).eq.'END') then
      print *, 'Early end of temp peptide file at line: ', k
      deallocate(peptide_pos) 
      close(unit=funit,iostat=ierr,status='keep')
      close(unit=funit2,iostat=ierr,status='keep')
      close(unit=funit3,iostat=ierr,status='keep')
      close(unit=funit4,iostat=ierr,status='keep')
      close(unit=funit5,iostat=ierr,status='keep')
      close(unit=funit6,iostat=ierr,status='keep')
      close(unit=funit7,iostat=ierr,status='keep')
      stop
    end if
    write(unit=funit4,fmt='(a)',iostat=ierr) &
  &    trim(word_string(13))//' '//trim(word_string(7))//' '// &
  &    trim(word_string(8))//' '//trim(word_string(9)) 
    write(unit=funit5,fmt='(a)',iostat=ierr) &
  &    trim(word_string(13))//' '//trim(word_string(7))//' '// &
  &    trim(word_string(8))//' '//trim(word_string(9)) 
    read(word_string(7),*) atom_pos(1)
    read(word_string(8),*) atom_pos(2)
    read(word_string(9),*) atom_pos(3)
    atom_pos_bohr = atom_pos/bohr
    call element_to_lower(word_string(13),str_element)
    write(unit=funit7,fmt='(f9.4,x,f9.4,x,f9.4,x,a2)',iostat=ierr) &
  &    atom_pos_bohr,str_element
    if (.not.selective) then
      read(word_string(7),*) peptide_pos(1,k)
      read(word_string(8),*) peptide_pos(2,k)
      read(word_string(9),*) peptide_pos(3,k)
      read(word_string(13),*) peptide_element(k)
    end if
  end do
  close(unit=funit2,iostat=ierr,status='delete')
  close(unit=funit4,iostat=ierr,status='keep')

!! Block to read in temp2 atoms data.
  open(unit=funit6,file=temp2_name,action='read',iostat=ierr)
  if (selective) then
    do l = 1, header_lines
      read(unit=funit6,fmt='(a)',iostat=ierr) file_string2
      if (ierr.ne.0) then
        print *, 'Problem reading temp2 file info, line: ', l
        deallocate(peptide_pos) 
        close(unit=funit,iostat=ierr,status='keep')
        close(unit=funit3,iostat=ierr,status='keep')
        close(unit=funit5,iostat=ierr,status='keep')
        close(unit=funit6,iostat=ierr,status='keep')
        close(unit=funit7,iostat=ierr,status='keep')
        stop
      end if
    end do
    do k = 1, n
      read(unit=funit6,fmt='(a)',iostat=ierr) file_string2
      if (ierr.ne.0) then
        print *, 'Problem reading temp2 peptide atom index: ', k
        deallocate(peptide_pos) 
        close(unit=funit,iostat=ierr,status='keep')
        close(unit=funit3,iostat=ierr,status='keep')
        close(unit=funit5,iostat=ierr,status='keep')
        close(unit=funit6,iostat=ierr,status='keep')
        close(unit=funit7,iostat=ierr,status='keep')
        stop
      end if
      call parse_PDB_line(file_string2,word_string,n_words)
      if (n_words.le.0) then
        print *, 'Problem reading temp2 peptide atom index: ', k
        print *, file_string2
        deallocate(peptide_pos) 
        close(unit=funit,iostat=ierr,status='keep')
        close(unit=funit3,iostat=ierr,status='keep')
        close(unit=funit5,iostat=ierr,status='keep')
        close(unit=funit6,iostat=ierr,status='keep')
        close(unit=funit7,iostat=ierr,status='keep')
        stop
      end if
      if (trim(word_string(1)).eq.'END') then
        print *, 'Early end of temp2 peptide file at line: ', k
        deallocate(peptide_pos) 
        close(unit=funit,iostat=ierr,status='keep')
        close(unit=funit3,iostat=ierr,status='keep')
        close(unit=funit5,iostat=ierr,status='keep')
        close(unit=funit6,iostat=ierr,status='keep')
        close(unit=funit7,iostat=ierr,status='keep')
        stop
      end if
      read(word_string(7),*) peptide_pos(1,k)
      read(word_string(8),*) peptide_pos(2,k)
      read(word_string(9),*) peptide_pos(3,k)
      read(word_string(13),*) peptide_element(k)
    end do
  end if
  close(unit=funit6,iostat=ierr,status='delete')

!! Exiting the peptide loop above leaves:
!!  the index j at the first non-peptide atom,
!!  the file_string contains the information about 
!!   this first non-peptide atom
!! the word_string has been destroyed and is regenerated below
!! This first atom needs to be processed before proceeding
!!  to read more data!
  call parse_PDB_line(file_string,word_string,n_words)
  if (n_words.le.0) then
    print *, 'Problem reading first non-peptide atom index: ',j
    print *, file_string
    deallocate(peptide_pos) 
    close(unit=funit,iostat=ierr,status='keep')
    close(unit=funit3,iostat=ierr,status='keep')
    close(unit=funit5,iostat=ierr,status='keep')
    close(unit=funit7,iostat=ierr,status='keep')
    stop
  end if

  n_shell_atoms = 0
  l = 0
  m = n_peptide_atoms
  do
!! word_string filled before loop or in previous loop
    read(word_string(7),*) atom_pos(1)
    read(word_string(8),*) atom_pos(2)
    read(word_string(9),*) atom_pos(3)
    atom_of_interest = .false.
    do k = 1, n
      if (peptide_element(k).eq.'O') then
        threshold = cutoff_O
      else if (peptide_element(k).eq.'N') then
        threshold = cutoff_N
      else if (peptide_element(k).eq.'S') then
        threshold = cutoff_S
      else
        threshold = cutoff
      end if
      distance_vec = atom_pos(:) - peptide_pos(:,k)
      distance = dot_product(distance_vec,distance_vec)
      distance = sqrt(distance)
      if ((threshold-distance).ge.eps_single) then
        atom_of_interest = .true.
        exit
      end if
    end do
    if (atom_of_interest) then
      m = m + 1
!! Special treatment of complex ions, needs to be hard coded
!!  for each ion to read and write in the entire complex.
      if (trim(word_string(4)).eq.'TIP3') then
!! WRITE THE FIRST ATOM (OXYGEN) OF THE COMPLEX
        write(unit=funit3,fmt='(a)',iostat=ierr) file_string
        write(unit=funit5,fmt='(a)',iostat=ierr) &
  &        trim(word_string(13))//' '//trim(word_string(7))//' '// &
  &        trim(word_string(8))//' '//trim(word_string(9)) 
        atom_pos_bohr = atom_pos/bohr
        write(unit=funit7,fmt='(f9.4,x,f9.4,x,f9.4,x,a2)',iostat=ierr) &
  &        atom_pos_bohr,'o '
        read(unit=funit,fmt='(a)',iostat=ierr) file_string
        if (ierr.ne.0) then
          print *, 'Unexpected problem reading TIP3 H1'
          print *, 'at index',j+1
          exit
        end if
        write(unit=funit3,fmt='(a)',iostat=ierr) file_string
        call parse_PDB_line(file_string,word_string,n_words)
        if (n_words.le.0) then
          print *, 'Unexpected problem parsing TIP3 H1'
          print *, 'at index',j+1
          exit
        end if
        write(unit=funit5,fmt='(a)',iostat=ierr) &
  &        trim(word_string(13))//' '//trim(word_string(7))//' '// &
  &        trim(word_string(8))//' '//trim(word_string(9)) 
        read(word_string(7),*) atom_pos(1)
        read(word_string(8),*) atom_pos(2)
        read(word_string(9),*) atom_pos(3)
        atom_pos_bohr = atom_pos/bohr
        write(unit=funit7,fmt='(f9.4,x,f9.4,x,f9.4,x,a2)',iostat=ierr) &
  &        atom_pos_bohr,'h '
        read(unit=funit,fmt='(a)',iostat=ierr) file_string
        if (ierr.ne.0) then
          print *, 'Unexpected problem reading TIP3 H2'
          print *, 'at index',j+2
          exit
        end if
        write(unit=funit3,fmt='(a)',iostat=ierr) file_string
        call parse_PDB_line(file_string,word_string,n_words)
        if (n_words.le.0) then
          print *, 'Unexpected problem parsing TIP3 H2'
          print *, 'at index',j+2
          exit
        end if
        write(unit=funit5,fmt='(a)',iostat=ierr) &
  &        trim(word_string(13))//' '//trim(word_string(7))//' '// &
  &        trim(word_string(8))//' '//trim(word_string(9)) 
        read(word_string(7),*) atom_pos(1)
        read(word_string(8),*) atom_pos(2)
        read(word_string(9),*) atom_pos(3)
        atom_pos_bohr = atom_pos/bohr
        write(unit=funit7,fmt='(f9.4,x,f9.4,x,f9.4,x,a2)',iostat=ierr) &
  &        atom_pos_bohr,'h '
        n_shell_atoms = n_shell_atoms + 3
        j = j + 2
        m = m + 2
      else
        write(unit=funit3,fmt='(a)',iostat=ierr) file_string
        l = l + 1
        print *, 'Non-water species ',trim(word_string(13)) ,' near peptide! Index:', j
        print *, trim(word_string(13))// &
  &       ' '//trim(word_string(7))//' '// &
  &       trim(word_string(8))//' '//trim(word_string(9)) 
      end if
    else
!! Special treatment of complex ions, needs to be hard coded
!!  to read and skip all atoms if first atom does not match
      if (trim(word_string(4)).eq.'TIP3') then
        read(unit=funit,fmt='(a)',iostat=ierr) file_string
        if (ierr.ne.0) then
          print *, 'Unexpected problem reading TIP3 H1'
          print *, 'at index',j+1
          exit
        end if
        read(unit=funit,fmt='(a)',iostat=ierr) file_string
        if (ierr.ne.0) then
          print *, 'Unexpected problem reading TIP3 H2'
          print *, 'at index',j+2
          exit
        end if
        j = j + 2
      end if
    end if
!!  move index to read next line
    j = j + 1
    read(unit=funit,fmt='(a)',iostat=ierr) file_string
    if (ierr.ne.0) then
      print *, 'Problem reading atom index: ',j
      print *, 'Or encountered end of input file.'
      exit
    end if
    call parse_PDB_line(file_string,word_string,n_words)
    if (n_words.le.0) then
      print *, 'Problem reading atom index: ',j
      print *, file_string
      exit
    end if
    if (trim(word_string(1)).eq.'END') then
      print *, 'Last line of file: ',j
      exit
    end if
    if (trim(word_string(1)).eq.'TER') then
      print *, 'Termination line of file: ',j
      exit
    end if
  end do

  print *, 'number of water solvent atoms selected: ',n_shell_atoms
  print *, 'number of non-water solvent atoms found: ',l
  print *, 'Total number atoms found:', m
  print *, 'Total number atoms selected:', (n_peptide_atoms+n_shell_atoms)

!! Done with file reading/writing
  deallocate(peptide_pos) 
  close(unit=funit,iostat=ierr,status='keep')
  write(unit=funit3,fmt='(a)',iostat=ierr) 'END'
  close(unit=funit3,iostat=ierr,status='keep')
  write(unit=funit5,fmt='(a)',iostat=ierr) ''
  write(unit=funit5,fmt='(a)',iostat=ierr) ''
  close(unit=funit5,iostat=ierr,status='keep')
  write(unit=funit7,fmt='(a)',iostat=ierr) '$end'
  close(unit=funit7,iostat=ierr,status='keep')


!--------------------------------------------------------------------
!--------------------------------------------------------------------
end program extract
!--------------------------------------------------------------------
!--------------------------------------------------------------------

!! subroutine to parse atom lines in .pdb files
subroutine parse_PDB_line(file_string,word_string,n_words)

use basekinds
  implicit none

  character(len=80), intent(in) :: file_string
  character(len=8), intent(out) :: word_string(13)
  integer(kind_integer), intent(out) :: n_words


  word_string(1) = file_string(1:4)
  word_string(2) = adjustl(file_string(5:11))
  word_string(3) = adjustl(file_string(12:16))
  word_string(4) = adjustl(file_string(17:21))
  word_string(5) = adjustl(file_string(22:23))
  word_string(6) = adjustl(file_string(24:26))
  word_string(7) = adjustl(file_string(31:38))
  word_string(8) = adjustl(file_string(39:46))
  word_string(9) = adjustl(file_string(47:54))
  word_string(10) = adjustl(file_string(55:60))
  word_string(11) = adjustl(file_string(61:66))
  word_string(12) = adjustl(file_string(67:73))
  word_string(13) = adjustl(file_string(74:78))
  n_words = 13

end subroutine parse_PDB_line

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

