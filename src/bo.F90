module bo_mod

  use base, only : rc2, copyptr
  use atoms
  use reaxff_param_mod

contains
!--------------------------------------------------------------------------------------------
SUBROUTINE BOCALC(nlayer, atype, pos)
!--------------------------------------------------------------------------------------------
integer,intent(in) :: nlayer
real(8),intent(in) :: atype(NBUFFER), pos(NBUFFER,3)

integer :: ti,tj,tk

!$omp master
call system_clock(ti,tk)
!$omp end master

!--- calculate BO prime 
CALL BOPRIM()  

!--- calculate full BO
CALL BOFULL() 

!$omp master
call system_clock(tj,tk)
it_timer(6)=it_timer(6)+(tj-ti)
!$omp end master

CONTAINS

!--------------------------------------------------------------------------------------------
SUBROUTINE BOPRIM()
!--------------------------------------------------------------------------------------------
! Calculates the BOp(0:3,i,j) and the deltap(i). 
!--------------------------------------------------------------------------------------------
implicit none

integer :: n,i,j, j1,i1
integer :: ity, jty, inxn
integer :: c1,c2,c3
real(8) :: dr(3), dr2, arg_BOpij(3)

!--- initialize deltap(1:,1) to -Val(atype(i))
!$omp do 
do i=1, copyptr(6)
   ity = nint(atype(i))
   deltap(i,1) = -Val(ity) 
enddo
!$omp end do

!$omp do
do i=1, copyptr(6)
   ity = nint(atype(i))

   do j1=1, nbrlist(i,0)

     j = nbrlist(i,j1)

     if(j<i) then
        jty  = nint(atype(j))
        inxn = inxn2(ity, jty)

        i1 = nbrindx(i,j1)

        dr(1:3) = pos(i,1:3) - pos(j,1:3)
        dr2= sum(dr(1:3)*dr(1:3))

        if(dr2 <= rc2(inxn)) then

          arg_BOpij(1) = cBOp1(inxn)*dr2**pbo2h(inxn)
          arg_BOpij(2) = cBOp3(inxn)*dr2**pbo4h(inxn)
          arg_BOpij(3) = cBOp5(inxn)*dr2**pbo6h(inxn)

          bo(1:3,i,j1) = switch(1:3,inxn)*dexp( arg_BOpij(1:3) )

!<kn> Small modification exists in sigma-bond prime, see reac.f line 4444. sigma-bond prime is multiplied by
!<kn> (1.d0 + 1.d-4) here.  Later in original reaxff code, sigma-bond prime is subtracted by 
!<kn>  0.01*vpar30 resulting in the following subtractions in abo(i1), bo(nbon), bos(nbon), bosi(nbon). 
!<kn> However, this modification is only applied to the energy calculation, not to the force calculation.
!<kn> The subtraction by <cutoff_vpar30> is done after the derivative calculations so that the variables in the
!<kn> force-calc routines use the original BOp values and the ones in the energy-calc routines are the modified value. 
          bo(1,i,j1) = ( 1.d0 + cutoff_vpar30 )*bo(1,i,j1)

!<kn> If the total <BOp> before the subtraction is greater than <cutoff_vpar30>, 
!<kn> get "final" <BOp> value and its derivatives.
          if(sum(bo(1:3,i,j1))>cutoff_vpar30 ) then 

             dln_BOp(1,i,j1)=switch(1,inxn)*pbo2(inxn)*arg_BOpij(1)
             dln_BOp(2,i,j1)=switch(2,inxn)*pbo4(inxn)*arg_BOpij(2)
             dln_BOp(3,i,j1)=switch(3,inxn)*pbo6(inxn)*arg_BOpij(3)
             dln_BOp(1:3,i,j1) = dln_BOp(1:3,i,j1)/dr2
             dln_BOp(1:3,j,i1) = dln_BOp(1:3,i,j1)

             dBOp(i,j1) = sum( bo(1:3,i,j1)*dln_BOp(1:3,i,j1) )
             dBOp(j,i1) = dBOp(i,j1)
          
!<kn> After the derivative calculations are done, do the subtraction described above 
!<kn> which results in the difference of bond-order  between the energy-calc and the force-calc.
             bo(1,i,j1) = bo(1,i,j1) - cutoff_vpar30

             bo(0,i,j1) = sum( bo(1:3,i,j1) ) 
             bo(0:3,j,i1) = bo(0:3,i,j1)

!$omp atomic
             deltap(i,1) = deltap(i,1) + bo(0,i,j1)
!$omp atomic
             deltap(j,1) = deltap(j,1) + bo(0,i,j1)
          else
             dBOp(i,j1) = 0.d0
             dBOp(j,i1) = 0.d0
             bo(0:3,i,j1) = 0.d0
             bo(0:3,j,i1) = 0.d0
          endif
        endif

      endif
   enddo
enddo
!$omp end do
 
END SUBROUTINE

!--------------------------------------------------------------------------------------------
SUBROUTINE BOFULL()
!--------------------------------------------------------------------------------------------
!  Subroutine calculates the Bond Order and its derivatives
!--------------------------------------------------------------------------------------------
implicit none

integer :: c1,c2,c3,n
integer :: i,i1,j,j1,ity,jty,inxn
real(8) :: exp1, exp2, exp12
real(8) :: fn1, fn2, fn3, fn4, fn5
real(8) :: fn23, fn45, fn145, fn1145
real(8) :: u1ij, u1ji
real(8) :: u1ij_inv2, u1ji_inv2
real(8) :: Cf1Aij, Cf1Bij
real(8) :: exp_delt22
real(8) :: Cf1ij, Cf1ji
real(8) :: pboc34
real(8) :: BOpij_2, u45ij, u45ji
real(8) :: exph_45ij, exph_45ji
real(8) :: Cf45ij, Cf45ji
real(8) :: fn45_inv
real(8) :: Cf1ij_div1, Cf1ji_div1
real(8) :: BOp0, BOpsqr

real(8) :: exppboc1i,exppboc2i,exppboc1j,exppboc2j   !<kn>

!$omp do
do i=1, copyptr(6)
   ity = nint(atype(i))
   deltap(i,2) = deltap(i,1) + Val(ity) - Valval(ity) ! update for Mo
enddo
!$omp end do

!$omp do 
do i=1, copyptr(6)
   ity = nint(atype(i))
  
   exppboc1i = exp( -vpar1*deltap(i,1) )  !<kn>
   exppboc2i = exp( -vpar2*deltap(i,1) )  !<kn>

   do j1=1, nbrlist(i,0)

      j=nbrlist(i,j1)

      if(j<i) then

        jty = nint(atype(j))

        exppboc1j = exp( -vpar1*deltap(j,1) )  !<kn>
        exppboc2j = exp( -vpar2*deltap(j,1) )  !<kn>

        i1=nbrindx(i,j1)

        inxn = inxn2(ity,jty)

        fn2 = exppboc1i + exppboc1j                                   !<kn>
        fn3 = ( -1.d0/vpar2 )*log( 0.5d0*(exppboc2i + exppboc2j) )    !<kn>

        fn23 = fn2 + fn3

!--- keep BOp(0) for later use 
        BOp0=bo(0,i,j1)

!--- check the over coordination flag <ovc>. If it's zero, fn1 remains as one.
!--- <ovc> is either 1.d0 or 0.d0 in the given parameter file.
        fn1 = 0.5d0*( ( Val(ity)+fn2 )/(Val(ity)+fn23 ) + (Val(jty)+fn2)/(Val(jty)+fn23) )
        if(ovc(inxn) < 1.d-3) fn1 = 1.d0

        BOpsqr = bo(0,i,j1)*bo(0,i,j1)
        fn4 = 1.d0/(1.d0 +  dexp(-pboc3(inxn) * (pboc4(inxn) * BOpsqr - deltap(i,2) ) + pboc5(inxn) ) )
        fn5 = 1.d0/(1.d0 +  dexp(-pboc3(inxn) * (pboc4(inxn) * BOpsqr - deltap(j,2) ) + pboc5(inxn) ) )

!--- Corresponding to <ovc>, <v13cor> is a flag to modify <fn4> and <fn5>.
!--- Currently (3/1/05) only Al-Al interaction satisfies this condition, which gives no-correction to 
!--- pi and double pi bond-order prime.
        if(v13cor(inxn)<1.d-3) then
          fn4 = 1.d0
          fn5 = 1.d0
        endif

        fn45   = fn4*fn5
        fn145  = fn1*fn45
        fn1145 = fn1*fn145

!--- New Bond-Order definition
        BO(0,i,j1) = bo(0,i,j1) * fn145
        BO(2,i,j1) = bo(2,i,j1) * fn1145
        BO(3,i,j1) = bo(3,i,j1) * fn1145
        if(BO(0,i,j1) < 1.d-10) BO(0,i,j1) = 0.d0
        if(BO(2,i,j1) < 1.d-10) BO(2,i,j1) = 0.d0
        if(BO(3,i,j1) < 1.d-10) BO(3,i,j1) = 0.d0

!--- new sigma BO definition.
        BO(1,i,j1) = BO(0,i,j1) - BO(2,i,j1) - BO(3,i,j1) 

        BO(0:3,j,i1) = BO(0:3,i,j1)
 
!--- CALCULATION OF DERIVATIVE OF BOND ORDER
!--- all following comes from Coding Methodology section:
!--- part 1:

        u1ij = Val(ity) + fn23
        u1ji = Val(jty) + fn23
 
!--- part 2:
        u1ij_inv2 = 1d0/(u1ij*u1ij)
        u1ji_inv2 = 1d0/(u1ji*u1ji)
  
        Cf1Aij = 0.5d0 * fn3 * (u1ij_inv2 + u1ji_inv2)
        Cf1Bij =-0.5d0 * ( (u1ij - fn3)*u1ij_inv2 + (u1ji - fn3)*u1ji_inv2 )

!--- part 3:
        exp_delt22 = exppboc2i + exppboc2j
        Cf1ij = ( -Cf1Aij*pboc1(inxn)*exppboc1i ) + (Cf1Bij*exppboc2i)/(exp_delt22)
        Cf1ji = ( -Cf1Aij*pboc1(inxn)*exppboc1j ) + (Cf1Bij*exppboc2j)/(exp_delt22)

!--- part 4:
        pboc34 = pboc3(inxn) * pboc4(inxn)  !consider array calcd outside
        BOpij_2 = BOpsqr

        u45ij = pboc5(inxn) + pboc3(inxn)*deltap(i,2) - pboc34*BOpij_2
        u45ji = pboc5(inxn) + pboc3(inxn)*deltap(j,2) - pboc34*BOpij_2

!--- part 5:
        exph_45ij = exp(u45ij)
        exph_45ji = exp(u45ji)
        exp1 = 1.d0/(1.d0 + exph_45ij)
        exp2 = 1.d0/(1.d0 + exph_45ji)
        exp12 = exp1*exp2

        Cf45ij = -exph_45ij*exp12*exp1
        Cf45ji = -exph_45ji*exp12*exp2

!--- if following conditions, <ovc> and <v13cor>, are satisfied, correction terms
!--- fn1 and/or fn4 & fn5 are 1.d0 and their derivatives are zero.  
        if(ovc(inxn) < 1.d-3) then 
           Cf1ij = 0.d0
           Cf1ji = 0.d0
        endif

        if(v13cor(inxn) < 1.d-3) then
           Cf45ij = 0.d0
           Cf45ji = 0.d0
        endif

!--- part 6:
        fn45_inv = 1d0/fn45

        Cf1ij_div1 = Cf1ij/fn1
        Cf1ji_div1 = Cf1ji/fn1

        A0(i,j1) = fn145
        A1(i,j1) = -2d0*pboc34*BOp0*(Cf45ij + Cf45ji)*fn45_inv
        A2(i,j1) = Cf1ij_div1 + (pboc3(inxn)*Cf45ij*fn45_inv)
        A3(i,j1) = A2(i,j1) + Cf1ij_div1

        A0(j,i1) = A0(i,j1)
        A1(j,i1) = A1(i,j1)
        A2(j,i1) = Cf1ji_div1 + (pboc3(inxn)*Cf45ji*fn45_inv)
        A3(j,i1) = A2(j,i1) + Cf1ji_div1     

      endif !if(i<j)

   enddo ! do j1=1,nbrlist(i,0) loop end
enddo
!$omp end do

!--- Calculate delta(i):
!--- Initialize detal(i)
!$omp do
do i=1, copyptr(6)
   ity = nint(atype(i))
   delta(i) = -Val(ity) + sum( BO(0,i,1:nbrlist(i,0)) )
enddo
!$omp end do

END SUBROUTINE

END SUBROUTINE BOCALC

end module
