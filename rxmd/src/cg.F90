!---------------------------------------------------------------------------------
module CG
use atoms
!---------------------------------------------------------------------------------
integer,parameter :: CG_MaxMinLoop = 500
integer,parameter :: CG_MaxLineMinLoop = 100
integer,parameter :: CG_MaxBracketLoop = 10

real(8) :: CG_EPS= 1d-16 ! a check to emit warning message

contains

!---------------------------------------------------------------------------------
subroutine ConjugateGradient(atype,pos)
!---------------------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER),pos(3,NBUFFER)
real(8) :: f(3,NBUFFER),v(3,NBUFFER),q(NBUFFER)

real(8) :: p(3,NBUFFER) ! search direction
real(8) :: gold(3,NBUFFER), gnew(3,NBUFFER) ! old and new gradients

integer :: cgLoop

if(myid==0)  print'(a50)', 'start structural optimization.'

do cgLoop = 0, CG_MaxMinLoop-1

   call LineMinimization(atype,pos,gnew,p)


enddo 

call MPI_Finalize(ierr)
stop 'successfully finished structural optimization. '

end subroutine ConjugateGradient

!---------------------------------------------------------------------------------
subroutine LineMinimization(atype,pos,f,p)
! input: atom type, initial coordinate, and search direction.
! output: final coordinate and gradient at the coordinate. 
!---------------------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER),pos(3,NBUFFER),f(3,NBUFFER),q(NBUFFER),vdummy(1,1)
real(8),intent(in) :: p(3,NBUFFER)
integer :: lineMinLoop
real(8) :: stepl

call BracketSearchRange(atype,pos,p,stepl)

do lineMinLoop = 0, CG_MaxLineMinLoop-1
   call QEq(atype, pos, q)
   call FORCE(atype, pos, f, q)

!--- migrate atoms after positions are updated
   call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atype, pos, vdummy, f, q)
enddo

end subroutine LineMinimization

!---------------------------------------------------------------------------------
function MoveAtomAndEvaluateEnergy(atype,pos,p,stepl) result(potentialEnergy)
!---------------------------------------------------------------------------------
real(8),intent(in) :: atype(NBUFFER),pos(3,NBUFFER),p(3,NBUFFER),stepl
real(8) :: potentialEnergy

real(8) :: atypeTmp(NBUFFER),posTmp(3,NBUFFER),fTmp(3,NBUFFER),qTmp(NBUFFER)
real(8) :: vdummy(1,1) ! <- not used

posTmp(1:3,1:NATOMS)=pos(1:3,1:NATOMS)+stepl*p(1:3,1:NATOMS)
atypeTmp(1:NATOMS)=atype(1:NATOMS)

call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atypeTmp, posTmp, vdummy, fTmp, qTmp)
call QEq(atypeTmp, posTmp, qTmp)
call FORCE(atypeTmp, posTmp, fTmp, qTmp)

potentialEnergy=sum(PE(1:13))

return 
end function MoveAtomAndEvaluateEnergy

!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine BracketSearchRange(atype,pos,p,stepl) 
! input: atom type, initial coordinate, and search direction.
! output: step length to bracket an energy minimum along the search direction.
!---------------------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER),pos(3,NBUFFER),p(3,NBUFFER),stepl
integer :: bracketingLoop

do bracketingLoop = 0, CG_MaxBracketLoop-1
   
enddo 

end subroutine BracketSearchRange

!---------------------------------------------------------------------------------
subroutine GoldenSectionSearch(left,right)
! input: left and right boundaries of search range. 
!---------------------------------------------------------------------------------
implicit none
real(8) :: left,right
real(8) :: left2,right2
real(8) :: phi = 1.d0/1.61803398875d0 ! inverse of golden ratio





end subroutine GoldenSectionSearch

!---------------------------------------------------------------------------------
subroutine PolynomialFitSearch()
!---------------------------------------------------------------------------------
implicit none


end subroutine PolynomialFitSearch

!---------------------------------------------------------------------------------
subroutine NormalizeVector(v, vnorm) 
!---------------------------------------------------------------------------------
implicit none
real(8) :: v(3,NBUFFER)
real(8) :: vsum, vnorm 
integer :: i

vsum=0.d0
do i=1, NATOMS
   vsum = vsum + sum(v(1:3,i)*v(1:3,i))
enddo

call MPI_ALLREDUCE(vsum, vnorm, 1, mpi_double_precision, mpi_sum,  mpi_comm_world, ierr)
vnorm=sqrt(vnorm)

if(abs(vnorm)<CG_EPS) then
   if(myid==0) print'(a,es25.15)', &
    'WARNING: Norm of vector was found too small in NormalizeVector(): vnorm = ', vnorm
endif

do i=1, NATOMS
   v(1:3,i)=v(1:3,i)/vnorm
enddo

return
end subroutine NormalizeVector

end module  
