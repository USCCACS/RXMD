!----------------------------------------------------------------------------------------------------------------------
module cg
!----------------------------------------------------------------------------------------------------------------------
real(8),allocatable :: cg_pos(:,:), cg_q(:), cg_atype(:)
real(8),allocatable :: h(:,:), g(:,:)

!--- max number of CG iteration
INTEGER, PARAMETER :: CGMAX=400 
REAL(8), PARAMETER :: EPS=1.0d-10

!--- mnbrak
REAL(8), PARAMETER :: GOLD=1.618034d0, GLIMIT=100.0d0, TINY=1.0d-20

!--- dbrent
INTEGER, PARAMETER :: ITMAX=100 

!--- dlinmin
REAL(8), PARAMETER :: TOL=1.d-4
REAL(8), PARAMETER :: CG_DX=1.d-3

end module

!---------------------------------------------------------------------------------
SUBROUTINE conjugate_gradient(atype, pos, v, f, q) !called frprmn in NR
use atoms; use cg
!---------------------------------------------------------------------------------
  IMPLICIT NONE 

real(8) :: atype(NBUFFER), q(NBUFFER)
real(8) :: pos(3,NBUFFER),v(3,NBUFFER),f(3,NBUFFER)

  INTEGER :: iter 
  REAL(8) :: fret 
  INTEGER :: its,i
  REAL(8) ::dgg,fp,gam,gg, fnrm,dum

  if(myid==0) then
     write(6,'(a15,e13.2)') 'INFO: Start CG ', ftol 
     write(6,'(a20)') 'INFO: reset velocity'
  endif
  v(:,:)=0.d0

  call cg_init()

  call QEq(4)
  call Force()

!--- initialize gradient and search direction
  g(:,:)=f(:,:)
  h(:,:)=g(:,:) 
  gg=sum(g(1:3,1:NATOMS)*g(1:3,1:NATOMS))
  dum = gg
  call MPI_ALLREDUCE(dum, gg, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

  fp=1d10
  do its=1, CGMAX

!--- save current number of iteration to return
     iter=its 

     call dlinmin(fret)

!--- print out force norm
     fnrm=0.d0
     do i=1,NATOMS
       fnrm=fnrm+sum(f(1:3,i)*f(1:3,i))
     enddo 
     dum = fnrm
     call MPI_ALLREDUCE (dum, fnrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
     fnrm=sqrt(fnrm)
     if(myid==0) write(6,'(a,i5,2d25.13)') 'its:  fret, fnrm', its, fret,log(fnrm)

!--- check convergence condition
!     if (2.0d0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN 
     if (abs((fret-fp)/fp) <= ftol) then 
        if(myid==0) print*,'Energy is converged'
        call OUTPUT_cg()
        call MPI_FINALIZE(ierr)
        stop
     endif
!--- save old energy to compare next step energy
     fp=fret 


!--- update gradient with obtained force in dlinmin
     g(:,:)=f(:,:)

!--- get conjugate direction coefficient
     dgg=sum(g(1:3,1:NATOMS)*g(1:3,1:NATOMS))
     dum = dgg
     call MPI_ALLREDUCE(dum, dgg, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
     if (gg == 0.0) then
        if(myid==0) print'(a)','Norm of gradient is zero'
        call MPI_FINALIZE(ierr)
        stop
     endif
     gam=dgg/gg 

!--- update search direction
     h(:,:)=g(:,:)+gam*h(:,:) 

!--- save old gradient norm
     gg = dgg 
  end do

if(myid==0) print*,'Energy is not converged'
call OUTPUT_cg()
call MPI_FINALIZE(ierr)
stop

return 

Contains 

!-------------------------------------------------------------------------------
SUBROUTINE dlinmin(fret) 
use atoms; use cg
!-------------------------------------------------------------------------------
  IMPLICIT NONE 
  REAL(8), INTENT(OUT) :: fret 
  REAL(8) :: aa,bb,fa,fb,fu,xmin,uu, Ghhmax,hhmax,hhmin 

  aa=0.0d0
  hhmax=maxval(h(1:3,1:NATOMS))
  hhmin=minval(h(1:3,1:NATOMS))
  hhmin=abs(hhmin)
  if(hhmax<hhmin) hhmax=hhmin
  call MPI_ALLREDUCE(hhmax, Ghhmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX,  MPI_COMM_WORLD, ierr)

  aa=0.0d0
  uu=CG_DX/Ghhmax

  call mnbrak(aa,uu,bb,fa,fu,fb) 
  fret=dbrent(aa,uu,bb,TOL,xmin) 

  pos=pos+xmin*h

  call xu2xs()
  call LINKEDLIST()
  call COPYATOMS(0)
  call xs2xu()

  call FORCE()

END SUBROUTINE dlinmin

!-------------------------------------------------------------------------------
SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc) 
use atoms; use cg
!-------------------------------------------------------------------------------
  IMPLICIT NONE 
  REAL(8), INTENT(INOUT) :: ax,bx 
  REAL(8), INTENT(OUT) :: cx,fa,fb,fc 
  REAL(8) :: fu,qq,r,u,ulim 
  real(8) :: dum

  call df1dim(ax,fa,dum)
  call df1dim(bx,fb,dum)

  if (fb > fa) then 
     call swap(ax,bx) 
     call swap(fa,fb) 
  end if
  cx=bx+GOLD*(bx-ax) 
  call df1dim(cx,fc,dum)
  do 
     if (fb<fc) RETURN 
     r=(bx-ax)*(fb-fc) 
     qq=(bx-cx)*(fb-fa) 
     u=bx-((bx-cx)*qq-(bx-ax)*r)/(2.0d0*sign(max(abs(qq-r),TINY),qq-r)) 
     ulim=bx+GLIMIT*(cx-bx) 
     if ((bx-u)*(u-cx) > 0.0d0) then
        call df1dim(u,fu,dum)
        if (fu <fc) then
           ax=bx 
           fa=fb 
           bx=u 
           fb=fu 
           RETURN 
        else if (fu > fb) then 
           cx=u 
           fc=fu 
           RETURN
        end if
        u=cx+GOLD*(cx-bx)
        call df1dim(u,fu,dum)
     else if ((cx-u)*(u-ulim) > 0.0d0) then
        call df1dim(u,fu,dum)
        if (fu <fc) then 
           bx=cx 
           cx=u 
           u=cx+GOLD*(cx-bx) 
           call df1dim(u,fu,dum)
           call shft(fb,fc,fu,fu) 
        end if
     else if ((u-ulim)*(ulim-cx) >= 0.0d0) then 
        u=ulim 
        call df1dim(u,fu,dum)
     else 
        u=cx+GOLD*(cx-bx) 
        call df1dim(u,fu,dum)
     end if
     call shft(ax,bx,cx,u) 
     call shft(fa,fb,fc,fu)
  end do

END SUBROUTINE mnbrak

SUBROUTINE shft(a,b,c,d) 
    REAL(8), INTENT(OUT) :: a 
    REAL(8), INTENT(INOUT) :: b,c 
    REAL(8), INTENT(IN) :: d 
    a=b 
    b=c 
    c=d 
END SUBROUTINE shft

SUBROUTINE swap(a,b) 
    REAL(8), INTENT(INOUT) :: a,b
    REAL(8) :: temp 
    temp = a 
    a=b
    b=temp
END SUBROUTINE swap

SUBROUTINE mov3(a,b,c,d,e,f) 
   REAL(8), INTENT(IN) :: d,e,f 
   REAL(8), INTENT(OUT) :: a,b,c 
   a=d 
   b=e 
   c=f 
END SUBROUTINE mov3

!-------------------------------------------------------------------------------
FUNCTION dbrent(ax,bx,cx,tol0,xmin) 
  use atoms; use cg 
!-------------------------------------------------------------------------------
  IMPLICIT NONE 
  REAL(8), INTENT(IN) :: ax,bx,cx,tol0 
  REAL(8), INTENT(OUT) :: xmin 
  REAL(8) ::dbrent 
  real(8),parameter :: ZEPS=1.d-10
  INTEGER :: iter 
  REAL(8) ::a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,& 
       u,u1,u2,vv,w,x,xm
  real(8) :: dum
  LOGICAL :: ok1,ok2

  a=min(ax,cx) 
  b=max(ax,cx) 
  vv=bx 
  w=vv 
  x=vv 
  e=0.0d0 
  call df1dim(x,fx,dx) 
  fv=fx 
  fw=fx 
  dv=dx 
  dw=dx 
  do iter=1, ITMAX 
     xm=0.5d0*(a+b) 
     tol1=tol0*abs(x) + ZEPS 
     tol2=2.0d0*tol1 
     if (abs(x-xm) <= (tol2-0.5d0*(b-a))) exit 
     if (abs(e) > tol1) then 
        d1=2.0d0*(b-a)
        d2=d1 
        if (dw /= dx) d1=(w-x)*dx/(dx-dw) 
        if (dv /= dx) d2=(vv-x)*dx/(dx-dv) 
        u1=x+d1 
        u2=x+d2 
        ok1=((a-u1)*(u1-b) > 0.0d0) .and. (dx*d1 <= 0.0d0) 
        ok2=((a-u2)*(u2-b) > 0.0d0) .and. (dx*d2 <= 0.0d0) 
        olde=e 
        e=d 
        if (ok1 .or. ok2) then 
           if(ok1 .and. ok2) then 
              d=merge(d1,d2, abs(d1) <abs(d2)) 
           else 
              d=merge(d1,d2,ok1) 
           end if
           if(abs(d) <=abs(0.5d0*olde)) then 
              u=x+d 
              if (u-a <tol2 .or. b-u < tol2) d=sign(tol1,xm-x) 
           else
              e=merge(a,b, dx >= 0.0d0)-x 
              d=0.5d0*e 
           end if
        else 
           e=merge(a,b, dx >= 0.0d0)-x 
           d=0.5d0*e 
        end if
     else 
        e=merge(a,b, dx >= 0.0d0)-x 
        d=0.5d0*e 
     end if

     if (abs(d) >= tol1) then 
        u=x+d 
        call df1dim(u,fu,dum) 
     else 
        u=x+sign(tol1,d) 
        call df1dim(u,fu,dum) 
        if (fu>fx) exit 
     end if

     call df1dim(u,dum,du) 
     if (fu <= fx) then 
        if (u>= x) then 
           a=x 
        else 
           b=x 
        end if
        call mov3(vv,fv,dv,w,fw,dw) 
        call mov3(w,fw,dw,x,fx,dx) 
        call mov3(x,fx,dx,u,fu,du) 
     else 
        if (u< x) then 
           a=u 
        else 
           b=u 
        end if
        if (fu <= fw.or. w== x) then 
           call mov3(vv,fv,dv,w,fw,dw) 
           call mov3(w,fw,dw,u,fu,du) 
        else if (fu <= fv .or. vv == x.or. vv== w) then 
           call mov3(vv,fv,dv,u,fu,du) 
        end if
     end if
  end do

  if (iter>ITMAX)  print*,'dbrent: exceeded maximum iterations'
  xmin=x 
  dbrent=fx 

END FUNCTION dbrent


!----------------------------------------------------------------------------------------------------------------------
subroutine df1dim(dx,fx,dfx)
use atoms; use cg
!----------------------------------------------------------------------------------------------------------------------
implicit none
real(8) :: dx,fx,dfx  
real(8) :: temp
integer :: i, NATOMS_old


!--- save the original position and atomtype
cg_pos(1:3,1:NATOMS)=pos(1:3,1:NATOMS)
cg_q(1:NATOMS)=q(1:NATOMS)
cg_atype(1:NATOMS)=atype(1:NATOMS)
NATOMS_old = NATOMS

!--- modify atom coord. with <dxh>
pos(1:3,1:NATOMS) = pos(1:3,1:NATOMS) + dx*h(1:3,1:NATOMS)

call xu2xs()
call LINKEDLIST()
call COPYATOMS(0)
call xs2xu()

!call QEq(4)
call FORCE()
PE(0)=sum(PE(1:13))
call MPI_ALLREDUCE (PE, GPE, size(PE), MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
GPE(:)=GPE(:)/GNATOMS

!--- get the estimate of gradient along <h> with the new atom coordinate
!--- <g> = <f> is gradient with negative sign 
dfx=0.d0
do i=1, NATOMS
   dfx = dfx - sum(f(1:3,i)*h(1:3,i))
enddo
!--- get the total gradient  
temp=dfx
call MPI_ALLREDUCE(temp, dfx, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

fx=GPE(0)

!--- Let <pos> and <atype> be the original value
NATOMS = NATOMS_old
pos(1:3,1:NATOMS)=cg_pos(1:3,1:NATOMS)
q(1:NATOMS)=cg_q(1:NATOMS)
atype(1:NATOMS)=cg_atype(1:NATOMS)

end subroutine
!----------------------------------------------------------------------------------------------------------------------
subroutine cg_init()
use atoms; use cg
!----------------------------------------------------------------------------------------------------------------------
allocate(cg_pos(3,NBUFFER),cg_q(NBUFFER),cg_atype(NBUFFER), stat=ast)
allocate(h(3,NBUFFER),g(3,NBUFFER), stat=ast)

end subroutine

!----------------------------------------------------------------------------------------
subroutine OUTPUT_cg()
use atoms
!----------------------------------------------------------------------------------------
implicit none
integer :: i
character(6) :: a6
character(9) :: a9

write(a6(1:6),'(i6.6)') myid
write(a9(1:9),'(i9.9)') nstep + current_step

open(10,file="DAT/rxff"//a6//"-"//a9//".bin",form="unformatted",access="stream")
write(10) NATOMS
write(10) current_step + ntime_step
write(10) lata,latb,latc,lalpha,lbeta,lgamma
do i=1,NATOMS
   write(10) pos(1:3,i)
   write(10) v(1:3,i)
   write(10) q(i)
   write(10) atype(i)
enddo
close(10)

end subroutine

END SUBROUTINE conjugate_gradient
