!---------------------------------------------------------------------------------
module CG
use atoms
!---------------------------------------------------------------------------------
integer,parameter :: CG_MaxMinLoop = 500
integer,parameter :: CG_MaxLineMinLoop = 100
integer,parameter :: CG_MaxBracketLoop = 20
integer,parameter :: CG_MaxGSLoop = 100


!--- Wolfe condition parameters
real(8),parameter :: CG_WC1 = 1d-4
real(8),parameter :: CG_WC2 = 0.1d0

!--- golden section tolerance
real(8),parameter :: CG_GStol = 1d-6

!--- conjugate gradient tolerance. ftol in rxmd.in
real(8) :: CG_tol = 1d-4

real(8),parameter :: CG_EPS= 1d-16 ! a check to emit warning message

contains

!---------------------------------------------------------------------------------
subroutine ConjugateGradient(atype,pos)
!---------------------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER),pos(NBUFFER,3)
real(8) :: f(NBUFFER,3),v(NBUFFER,3),q(NBUFFER)

real(8) :: p(NBUFFER,3) ! search direction
real(8) :: gold(NBUFFER,3),gnew(NBUFFER,3) ! old and new gradients
real(8) :: GPE(0:13),GPEold,GPEnew

real(8) :: vnorm, stepl, beta1,beta2,beta3

integer :: cgLoop, i

CG_tol = ftol
v(:,:)=0.d0

if(myid==0) print'(a40,1x,es9.2)', NEW_LINE('A')//'Start structural optimization.', ftol

call QEq(atype, pos, q)
call FORCE(atype, pos, gnew, q)

!--- initialize search direction with gradient
p(1:3,1:NATOMS)=gnew(1:3,1:NATOMS)

PE(0)=sum(PE(1:13))
call MPI_ALLREDUCE(PE, GPE, size(PE), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
GPEnew=GPE(0)

!--- if no bracket range was found here, you are at the energy minimum already. 
call BracketSearchRange(atype,pos,p,stepl)

do cgLoop = 0, CG_MaxMinLoop-1

   call LineMinimization(atype,pos,p,gnew,stepl)
   gold(1:3,1:NATOMS)=gnew(1:3,1:NATOMS)

   call QEq(atype, pos, q)
   call FORCE(atype, pos, gnew, q)

   call OUTPUT(atype, pos, v, q, GetFileNameBase(cgLoop))

   GPEold=GPEnew
   PE(0)=sum(PE(1:13))
   call MPI_ALLREDUCE(PE, GPE, size(PE), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   GPEnew=GPE(0)

   if(abs(GPEnew-GPEold)<=CG_tol*GNATOMS) then
      if(myid==0) print'(a30,i6)','Energy converged.', cgLoop

      call OUTPUT(atype, pos, v, q, GetFileNameBase(cgLoop))
      exit
   endif

   beta1=DotProductVec3D(gold,gold,NATOMS)
   beta2=DotProductVec3D(gnew,gnew,NATOMS)
   beta3=DotProductVec3D(gnew,gold,NATOMS)

   if(myid==0) print'(a30,i6,4es15.5,2es20.10)', &
     'b1,b2,b3,(b2-b3)/b1: ',cgLoop,beta1,beta2,beta3,(beta2-beta3)/beta1,GPEnew,GPEold

   p(1:3,1:NATOMS) = (beta2-beta3)/beta1*p(1:3,1:NATOMS) + gnew(1:3,1:NATOMS)

   call BracketSearchRange(atype,pos,p,stepl)

enddo 

call MPI_Finalize(ierr)
stop 'successfully finished structural optimization. '

end subroutine ConjugateGradient

!---------------------------------------------------------------------------------
subroutine BracketSearchRange(atype,pos,p,stepl) 
! input: atom type, initial coordinate, and search direction.
! output: step length to bracket an energy minimum along the search direction.
!---------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: atype(NBUFFER),pos(NBUFFER,3),p(NBUFFER,3)
real(8),intent(out) :: stepl
real(8) :: vdummy(NBUFFER,3), qdummy(NBUFFER)
integer :: bracketingLoop
logical :: Elower, WolfeC1, WolfeC2

real(8) :: PE

if(myid==0) print'(a40)', NEW_LINE('A')//'Start BracketSearchRange()'

stepl=1d-2/GNATOMS; Elower=.true.; WolfeC1=.true.; WolfeC2=.true.

do bracketingLoop = 0, CG_MaxBracketLoop-1

   stepl = stepl*2

   PE = EvaluateEnergyWithStep(atype,pos,p,stepl)
   call WolfeConditions(atype,pos,p,stepl,Elower,WolfeC1,WolfeC2)

   if(myid==0) print'(a30,es15.5,es25.15, 3l3)', &
      'stepl,PE,Elow,Wolfe1,Wolfe2: ', stepl, PE, Elower, WolfeC1, WolfeC2

   if(.not.WolfeC1 .or. .not.WolfeC1) then
      if(myid==0) print'(a30,es15.5,a1)', 'bracket has been found: ', stepl
      return 
   endif

enddo 

if(myid==0) print'(a)', &
'bracket was not found. saving the last configuration and terminating structural optimization'
call OUTPUT(atype, pos, vdummy, qdummy, "nobracket")

stop

end subroutine BracketSearchRange

!---------------------------------------------------------------------------------
subroutine WolfeConditions(atype,pos,p,stepl,isLowerEnergy,isArmijoRule,isCurvature)
! input: atom type, position, search direction and step length
! output: two bools for the Wolfe conditions
! TODO: This function is very similar to EvaluateEnergyWithStep because the function 
! doesn't return force and new NATOM. Can be simplified. 
!---------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: atype(NBUFFER),pos(NBUFFER,3),p(NBUFFER,3),stepl
logical,intent(out) :: isLowerEnergy,isArmijoRule,isCurvature

real(8) :: atypeTmp(NBUFFER),posTmp(NBUFFER,3),pTmp(NBUFFER,3),qTmp(NBUFFER)
! FIXME: here we don't really need v but COPYATOM(MOVE) requires it for the move mode. 
!        thus, I am allocating a 3xNBUFFER dummy array here. a better implementation needed. 
real(8) :: vdummy(NBUFFER,3) 
real(8) :: GPE(0:13), GPEbefore, GPEafter, fbefore(NBUFFER,3), fafter(NBUFFER,3)
real(8) :: pDotdF, pDotdFShift
integer :: NATOMSTmp

! Evaluate df(x) and f(x)
call QEq(atype, pos, qTmp)
call FORCE(atype, pos, fbefore, qTmp)

PE(0)=sum(PE(1:13))
call MPI_ALLREDUCE(PE, GPE, size(PE), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
GPEbefore=GPE(0)

! Evaluate df(x+alpha*p) and f(x+alpha*p)
posTmp(1:NATOMS,1:3)=pos(1:NATOMS,1:3)+stepl*p(1:NATOMS,1:3)
atypeTmp(1:NATOMS)=atype(1:NATOMS)

! FIXME local number of atoms will change after COPYATOM(MOVE). 
! Need to retrieve the origianl value consistent with the coordinates before the move.
NATOMSTmp = NATOMS

call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atypeTmp, posTmp, vdummy, fafter, qTmp)
call QEq(atypeTmp, posTmp, qTmp)
call FORCE(atypeTmp, posTmp, fafter, qTmp)

PE(0)=sum(PE(1:13))
call MPI_ALLREDUCE(PE, GPE, size(PE), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
GPEafter=GPE(0)

isLowerEnergy = GPEafter < GPEbefore

pDotdF = DotProductVec3D(p,fbefore,NATOMS)
isArmijoRule = GPEafter <= GPEbefore + pDotdF*CG_WC1*stepl

! the local number of atoms of f'(x+a) and p may differ after the shift, 
! cannot take their vector dot-product without moving the search vector. 

NATOMS = NATOMSTmp 
posTmp(1:NATOMS,1:3)=pos(1:NATOMS,1:3)+stepl*p(1:NATOMS,1:3)
pTmp(1:NATOMS,1:3)=p(1:NATOMS,1:3)
atypeTmp(1:NATOMS)=atype(1:NATOMS)
call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0],atypeTmp,posTmp,pTmp,vdummy,qTmp)

pDotdFShift = DotProductVec3D(pTmp,fafter,NATOMS)
isCurvature = pDotdFShift >= CG_WC2*pDotdF

! recover the original NATOM
NATOMS = NATOMSTmp 

end subroutine WolfeConditions

!---------------------------------------------------------------------------------
subroutine LineMinimization(atype,pos,p,g,stepl)
! Here, we perform line minimization based on golden section search. After obtaining,
! step length, update atom positions and migrate them if they move out of MD box. 
! We also migrate gradient and search vectors according to the position migration 
! in order to perform dot product of old and new gradient vectors. 
! input: atom type, initial coordinate, and search direction.
! output: updated coordinate and associated gradient vector. 
!---------------------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER),pos(NBUFFER,3),g(NBUFFER,3),p(NBUFFER,3),q(NBUFFER)
real(8) :: atypeTmp(NBUFFER),posTmp(NBUFFER,3),fdummy(1,1)
integer :: lineMinLoop, NATOMSTmp
real(8) :: stepl, step0

if(myid==0) print'(a40)', NEW_LINE('A')//'Start LineMinimization()'

step0=0.d0

call GoldenSectionSearch(atype,pos,p,step0,stepl)

! First, migrate gradient vector g.
NATOMStmp = MigrateVec3D(pos,p,g,stepl)

! Then, migrate atom type, position, and search vector; atype, pos, p.
pos(1:NATOMS,1:3)=pos(1:NATOMS,1:3)+stepl*p(1:NATOMS,1:3)
call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atype, pos, p, fdummy, q)

return
end subroutine LineMinimization

!---------------------------------------------------------------------------------
subroutine GoldenSectionSearch(atype,pos,p,ax,dx)
! ax,dx: left and right boundaries of search range. 
!---------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: atype(NBUFFER),pos(NBUFFER,3),p(NBUFFER,3)
real(8) :: ax,bx,cx,dx,PEbx,PEcx
real(8) :: ratio = 1.d0/1.61803398875d0 ! inverse of golden ratio

integer :: GSLoop

if(myid==0) print'(a30)', 'start golden section step.'

bx=dx-(dx-ax)*ratio
cx=ax+(dx-ax)*ratio
PEbx=EvaluateEnergyWithStep(atype,pos,p,bx)
PEcx=EvaluateEnergyWithStep(atype,pos,p,cx)

do GSLoop = 0, CG_MaxLineMinLoop-1

   if(myid==0) print'(a30,4es15.5,1x,2es25.15)', &
      'ax,bx,cx,dx,PEbx,PEcx: ', ax,bx,cx,dx,PEbx,PEcx

   if(abs(ax-dx)<=CG_GStol/GNATOMS) then
      if(myid==0) print'(a30)', 'golden section step finished.'
      exit
   endif

   if(PEbx<PEcx) then
      dx=cx
   else
      ax=bx
   endif

   bx=dx-(dx-ax)*ratio
   cx=ax+(dx-ax)*ratio
   PEbx=EvaluateEnergyWithStep(atype,pos,p,bx)
   PEcx=EvaluateEnergyWithStep(atype,pos,p,cx)
enddo

end subroutine GoldenSectionSearch

!---------------------------------------------------------------------------------
subroutine PolynomialFitSearch()
!To be implemented.
!---------------------------------------------------------------------------------
implicit none

end subroutine PolynomialFitSearch

!---------------------------------------------------------------------------------
function MigrateVec3D(pos, vec, dir, stepl) result(newNATOMS)
!TODO: come up a better way to migrate vectors
!---------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: pos(NBUFFER,3),stepl,dir(NBUFFER,3)
real(8) :: vec(NBUFFER,3)
real(8) :: atypedummy(NBUFFER),posTmp(NBUFFER,3),fdummy(1,1),qdummy(NBUFFER)
integer :: NATOMSTmp, newNATOMS

!--- keep current NATOMS
NATOMSTmp=NATOMS

posTmp(1:NATOMS,1:3)=pos(1:NATOMS,1:3)+stepl*dir(1:NATOMS,1:3)
call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atypedummy, posTmp, vec, fdummy, qdummy)

!-- this NATOMS is consistent with the migrated vector.
newNATOMS=NATOMS

!--- save the original NATOMS that is consistent with input pos. 
NATOMS=NATOMSTmp

return
end function MigrateVec3D


!---------------------------------------------------------------------------------
function DotProductVec3D(v1, v2, Nelems) result (ret)
!---------------------------------------------------------------------------------
implicit none
integer,intent(in) :: Nelems
real(8),intent(in) :: v1(Nelems,3), v2(Nelems,3) 
integer :: i
real(8) :: vsum, ret

vsum = 0.d0
do i=1, Nelems
   vsum = vsum + sum(v1(i,1:3)*v2(i,1:3))
enddo

call MPI_ALLREDUCE(vsum, ret, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

return
end function DotProductVec3D


!---------------------------------------------------------------------------------
subroutine NormalizeVec3D(v, vnorm) 
!---------------------------------------------------------------------------------
implicit none
real(8) :: v(NBUFFER,3)
real(8) :: vsum, vnorm 
integer :: i

vnorm=sqrt(DotProductVec3D(v,v,NATOMS))

if(abs(vnorm)<CG_EPS) then
   if(myid==0) print'(a,es25.15)', &
    'WARNING: Norm of vector was found too small in NormalizeVector(): vnorm = ', vnorm
endif

v(1:NATOMS,1:3)=v(1:NATOMS,1:3)/vnorm

return
end subroutine NormalizeVec3D

!---------------------------------------------------------------------------------
function EvaluateEnergyWithStep(atype,pos,p,stepl) result(potentialEnergy)
!---------------------------------------------------------------------------------
real(8),intent(in) :: atype(NBUFFER),pos(NBUFFER,3),p(NBUFFER,3),stepl
real(8) :: potentialEnergy

real(8) :: atypeTmp(NBUFFER),posTmp(NBUFFER,3),fTmp(NBUFFER,3),qTmp(NBUFFER)
! TODO: v is dummy but COPYATOM(MOVE) moves it. Thus needs to allocate 3xNBUFFER array here.
real(8) :: vdummy(NBUFFER,3) 
real(8) :: GPE(0:13)
integer :: NATOMSTmp

posTmp(1:NATOMS,1:3)=pos(1:NATOMS,1:3)+stepl*p(1:NATOMS,1:3)
atypeTmp(1:NATOMS)=atype(1:NATOMS)

! TODO: local number of atoms will change after COPYATOM(MOVE). 
! Need to retrieve the origianl value consistent with the coordinates before the move.
NATOMSTmp = NATOMS

call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atypeTmp, posTmp, vdummy, fTmp, qTmp)
call QEq(atypeTmp, posTmp, qTmp)
call FORCE(atypeTmp, posTmp, fTmp, qTmp)

NATOMS = NATOMSTmp 

PE(0)=sum(PE(1:13))
call MPI_ALLREDUCE(PE, GPE, size(PE), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
potentialEnergy=GPE(0)

end function EvaluateEnergyWithStep

!---------------------------------------------------------------------------------


end module  
