    module m_random
	implicit none
	character(len=64), parameter, private :: fic_seed="my_seed.dat"
!
	contains
!
	subroutine reset_seed(iseed)
	implicit none
!	integer, intent(in) :: iseed
	integer :: iseed
	call random_seed(iseed)
	end subroutine reset_seed
!
        subroutine save_seed()
        implicit none
        integer :: n, I_unit_seed
        integer, dimension(:), allocatable :: last_seed
        call random_seed(size=n)
        allocate(last_seed(n))
        call random_seed(get=last_seed)
!!$ write(6,*) "last_seed=",last_seed
        call get_unit(I_unit_seed)

open(unit=I_unit_seed,file=trim(adjustl(fic_seed)),status="unknown",form="unformatted")
write(I_unit_seed) n
write(I_unit_seed) last_seed
close(I_unit_seed)
deallocate(last_seed)
end subroutine save_seed
!
subroutine load_seed()
implicit none
integer :: n, I_unit_seed
integer :: iseed
logical :: L_present
integer, dimension(:), allocatable :: last_seed
!
call get_unit(I_unit_seed)

open(unit=I_unit_seed,file=trim(adjustl(fic_seed)),status="old",form="unformatted",err=100)
!write(6,*) "Loading seed file"
read(I_unit_seed) n
allocate(last_seed(n))
read(I_unit_seed) last_seed
close(I_unit_seed)
call random_seed(put=last_seed)
deallocate(last_seed)
return
!
100 continue
write(6,*) "Creating seed"
iseed=0
call reset_seed(iseed)
end subroutine load_seed
!
subroutine get_unit(Num_Fich)
implicit none
integer :: Num_Fich
logical :: ouvert
ouvert=.true.
Num_Fich=10
do while (Num_Fich < 100)
Num_Fich=Num_Fich+1
inquire(unit=Num_Fich,opened=ouvert)
if (.not.ouvert) exit
enddo
if (ouvert) then
write(6,*) "Pas d unite logique libre"
stop 'get_unit'
endif
end subroutine get_unit
!
end module m_random
!
SUBROUTINE RANDOMIZE(HARVEST)
use m_random
IMPLICIT NONE
integer :: i


INTEGER :: SEED
REAL :: HARVEST
REAL, DIMENSION(4,4) :: HARVEYS

!
call load_seed()
!
CALL RANDOM_NUMBER(HARVEST)
!CALL RANDOM_NUMBER(HARVEYS)

1 FORMAT(4F6.2)
!DO i=1,4
!WRITE(6,1) harveys(i,1), harveys(i,2), harveys(i,3), harveys(i,4)
!END DO
!WRITE(6,1) harvest
!
call save_seed()
!
END SUBROUTINE RANDOMIZE