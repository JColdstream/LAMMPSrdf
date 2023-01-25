program rdf
use omp_lib

implicit none
include 'comrdf.inc'
integer(8):: iter

frame = 1
call start
open(unit = 110, file='test.xyz', status='unknown')

do while (lastframe .eq. 0 .or. frame .le. lastframe)
  if (mod(frame, 50) .eq. 0) write(*,*) 'FRAME : ', frame
  call readheader
  if (trajstatus .ne. 0) exit
  call readcoordinates
  call initialcom
  call iteratecom
  write(110, *) nanalyse+1
  write(110, *)
  write(110, *) 'Zn', com(:)
  do iter = 1, nanalyse
    write(110, *) 'C', atom(:, iter)
  enddo
  call radialdensityfunction
  frame = frame+1
enddo

! test .xyz file to check molecule rebuilding
! open(unit = 110, file='test.xyz', status='unknown')
! do iter = 1, nanalyse
!   write(110, *) 'C', atom(:, iter)
! enddo

write(*, *) 'Final frame: ', frame-1
call rdf_output

contains

subroutine start

call get_command_argument(1, inputfile)
if (len_trim(inputfile) == 0) then
  write(*,*) 'Provide an input file.'
  call exit
endif

call readinput

! Open lammps trajectory file.
open(11, file = trajfile, status='old')
! get ntotal
call readheader

! Allocate arrays.
call readmasses
call nanalyse_setup
allocate(atom(3, nanalyse))
allocate(totalatom(3, ntotal))
allocate(atomtype(nanalyse))
allocate(totalatomtype(ntotal))
! Go back to the start of the file.
rewind(11)

! call skipframe
frame = frame+nskip

end subroutine start

! Reads analysis.input input file.
subroutine readinput
open(10, file = inputfile, status='old')
read(10, *)
read(10, *)
read(10, '(a)') trajfile
read(10, *) massfile
read(10, *) nskip
read(10, *) lastframe
read(10, *) ntypes
read(10, '(a)') strtypesanalyse
read(10, *) binwidth
typesanalyse = string_to_integers(strtypesanalyse, ",")
write(*,*) 'Trajectory file : ', trajfile


end subroutine readinput

! Reads the masses from the file given in the input.
subroutine readmasses
  integer(sp):: i, j

  allocate(mass(ntypes))
  if (massfile .eq. '0') then
    do i = 1, ntypes
      mass(i) = 1.0_dp
    enddo
  else
    open(12, file = massfile, status='old')
    read(12, *)
    read(12, *)
    do i = 1, ntypes
      read(12, *) j, mass(j)
    enddo
  endif
end subroutine readmasses

! Skips frames that don't need to be analysed.
subroutine skipframe
integer:: i
do i = 1, nskip*(ntotal+9)
  read(11, *)
enddo
write(*,*) "SKIPPED FRAMES : ", nskip
write(*,*)
end subroutine skipframe

! Reads in the box length, timestep, and natoms from the LAMMPS trajectory headers.
subroutine readheader
read(11, *, iostat = trajstatus) headertext(1)   
read(11, *, iostat = trajstatus) timestep
read(11, *, iostat = trajstatus) headertext(2)
read(11, *, iostat = trajstatus) ntotal
read(11, *, iostat = trajstatus) headertext(3)
read(11, *, iostat = trajstatus) xmin, xmax
read(11, *, iostat = trajstatus) ymin, ymax
read(11, *, iostat = trajstatus) zmin, zmax
read(11, *, iostat = trajstatus) headertext(4)

!if(natom  .ne. ntotal) stop 'mismatched number of atoms'

lx = xmax-xmin
ly = ymax-ymin
lz = zmax-zmin

! Edit pbc conditions if it needs to be used for a non-cubic box.
if ((lx .ne. ly) .or. (ly .ne. lz)) then
  write(*, *) 'Script is not suitable for non-cubic boxes.'
  call exit
endif

!write(*,*) timestep
!write(*,*) xmin, xmax
!write(*,*) ymin, ymax
!write(*,*) zmin, zmax

end subroutine readheader

subroutine nanalyse_setup
  integer(sp):: i, atom_type, nindex
  nanalyse = 0
  masstotal = 0.0_dp
  do i = 1, ntotal
    read(11, *) nindex, atom_type
    if (atom_type .gt. ntypes) then
      write(*,*) "Number of atom types mismatch." 
      call exit
    endif
    if ( any(typesanalyse .eq. atom_type) ) then
      masstotal = masstotal+mass(atom_type)
      nanalyse = nanalyse+1
    endif
  enddo
  ! Array to keep track of which atoms we need to calculate.
  allocate(molanalyse(nanalyse))
  nbins = (int((xmax-xmin)/binwidth)+1)
  allocate(rdist(ntypes, nbins))
end subroutine nanalyse_setup

subroutine readcoordinates
  integer(sp):: i, j, nindex, temp_atom_type
  real(dp):: tempx, tempy, tempz
  j = 1
  do i = 1, ntotal
    read(11, *, iostat = trajstatus)  nindex, temp_atom_type, tempx, tempy, tempz

    if (any(typesanalyse .eq. temp_atom_type)) then
      atom(:, j) = (/tempx, tempy, tempz/)
      atomtype(j) = temp_atom_type
      j = j+1
    endif
    totalatomtype(i) = temp_atom_type
    totalatom(:, i) = (/tempx, tempy, tempz/)
  enddo

 ! Zero coordinates.
 atom(1, :) = atom(1, :)-xmin
 atom(2, :) = atom(2, :)-ymin
 atom(3, :) = atom(3, :)-zmin
end subroutine readcoordinates

subroutine initialcom
  integer(8):: mol, i

  com = 0.0_dp
  do mol = 1, nanalyse
    do i = 1, 3
      com(i) = com(i) + atom(i, mol)*mass(atomtype(mol))
    enddo
  enddo
  com = com/masstotal
end subroutine initialcom

subroutine iteratecom
  integer(8):: i, mol
  real(dp), dimension(3):: temp_com, dxyz, com_diff
  real(dp):: com_diff_mag, drcom 
  logical:: final_com

  final_com = .False.

  do while (final_com .eqv. .False.)
    temp_com = com
    com = 0.0_dp
    do mol = 1, nanalyse
      do i = 1, 3
        ! dxyz(i) = dxyz(i) + (lx*anint(dxyz(i)/lx))
        com(i) = com(i) + (atom(i, mol) - lx*anint((atom(i, mol)-temp_com(i))/lx)) * mass(atomtype(mol))
        ! com(i) = com(i) + atom(i, mol) - lx*anint((atom(i, mol)-temp_com(i))/lx)
      enddo
    enddo
    com = com/masstotal
    !com = com/nanalyse
    do i = 1, 3
      if (com(i) .gt. lx) then
        com(i) = com(i) - lx
      elseif (com(i) .lt. 0.0) then
        com(i) = com(i) + lx
      endif
    enddo
    ! write(*,*) com
    com_diff = com-temp_com
    ! check how close the new centre of geometry is to the previous iteration
    com_diff_mag = 0.0_dp
    do i = 1, 3
      com_diff_mag = com_diff_mag+com_diff(i)**2
    enddo
    if (com_diff_mag .eq. 0.0_dp) then
      final_com = .True.
    endif
  enddo

  do mol = 1, nanalyse
    do i = 1, 3
      drcom = atom(i, mol) - com(i)
      if ( abs(drcom) .gt. lx/2 ) then
        atom(i, mol) = atom(i, mol) - lx*anint(drcom/lx)
      endif
    enddo
  enddo

  do mol = 1, ntotal
    do i = 1, 3
      drcom = totalatom(i, mol) - com(i)
      if ( abs(drcom) .gt. lx/2 ) then
        totalatom(i, mol) = totalatom(i, mol) - lx*anint(drcom/lx)
      endif
    enddo
  enddo

end subroutine iteratecom

subroutine radialdensityfunction
  integer(8)::  i, mol, moli, molj, rbin
  real(dp):: drsq, dr, drcom, vshell

  !$OMP PARALLEL DO &
  !$OMP SCHEDULE(DYNAMIC) &
  !$OMP DEFAULT(NONE) &
  !$OMP SHARED(nq, q, lx, atom, nanalyse, b, totalatomtype, rg) &
  !$OMP PRIVATE(iq, moli, molj, drsq, dr, qrij, dxyz, i, tempq) &
  !$OMP REDUCTION(+:p)

  do mol = 1, ntotal
    drsq = 0.0_dp
    do i = 1, 3
      drsq = drsq + (totalatom(i, mol) - com(i))**2
    end do  
    rbin = int(sqrt(drsq)/binwidth+1)
    rdist(totalatomtype(mol), rbin) = rdist(totalatomtype(mol), rbin) + 1.0
  end do  
  !$OMP END PARALLEL DO

  rdistnrm = rdistnrm+1.0

end subroutine radialdensityfunction

subroutine rdf_output
  integer(8):: i, j, t
  real(dp):: vshell

  ! calculate the density
  do t = 1, ntypes
    do i = 1, nbins
      vshell = (4.0/3.0) * pi * ((i*binwidth)**3 - ((i-1)*binwidth)**3)
      rdist(t, i) = rdist(t, i)/vshell 
    enddo
  enddo 

  open(unit = ndensf, file='ndensity.dat', status='unknown')
  open(unit = mdensf, file='massdensity.dat', status='unknown')
  write(ndensf, *) '#  r', ('    ', j, j = 1, ntypes)
  write(mdensf, *) '#  r', ('    ', j, j = 1, ntypes)

  do i = 1, nbins
    write(ndensf, '(f10.3, 999(e15.6))') (binwidth*(dble(i)-0.5)), &
      (rdist(j, i)/rdistnrm, j = 1, ntypes)
    write(mdensf, '(f10.3, 999(e15.6))') (binwidth*(dble(i)-0.5)), & 
      (rdist(j, i)*mass(j)/rdistnrm, j = 1, ntypes)
  enddo
  write(*,*) 'Number densities written to ndensity.dat'
  write(*,*) 'Number densities written to massdensity.dat'
  
end subroutine rdf_output

! thanks Rui!
function string_to_integers(str, sep) result(a)
    integer, allocatable:: a(:)
    character(*):: str
    character:: sep
    integer:: i, n_sep

    n_sep = 0
    do i = 2, len_trim(str)
      if (str(i:i)==sep .and. str(i-1:i-1)/=sep) then
        n_sep = n_sep+1
        str(i:i) = ','
       end if
    end do
    allocate(a(n_sep+1))
    read(str, *) a
end function

end program rdf
