!----------------------------------------------------------------------------------------------------------------------
subroutine FORCE(NBUFFER, atype, pos, f, q)
use parameters; use atoms 
!----------------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(in) :: NBUFFER
real(8) :: atype(NBUFFER), q(NBUFFER)
real(8) :: pos(3,NBUFFER), f(3,NBUFFER), vdummy(3,NBUFFER)

integer :: i, j, j1, ity, l(3),k
real(8) :: fsum(3), ss, vv(3), rr(3)
integer :: l2g
real(8) :: gasdev
real(8) :: dr(3)

ccbnd(:) = 0.d0
f(:,:) = 0.d0
fnb(:,:) = 0.d0

#ifdef STRESS
!--- stress components have to be transfered back to the original atoms, as the force components. 
astr(:,:) = 0.d0
#endif

!--- unscaled to scaled coordinate
call xu2xs(NBUFFER, pos)
call system_clock(i,k)
dr(1:3)=NMINCELL*lcsize(1:3)
call COPYATOMS(MODE_COPY,dr,NBUFFER,atype,pos,vdummy,f,q) 
call system_clock(j,k)
it_timer(4)=it_timer(4)+(j-i)

call system_clock(i,k)
call LINKEDLIST(NBUFFER, atype, pos)
call NBLINKEDLIST(NBUFFER, atype, pos)
call system_clock(j,k)
it_timer(3)=it_timer(3)+(j-i)
call xs2xu(NBUFFER, pos)

!--- scaled to unscaled coordinate

!$omp parallel default(shared), private(i,j,k)
!$omp sections

!$omp section

call system_clock(i,k)
call GetNonbondingPairList(NBUFFER, atype, pos)
call system_clock(j,k)
it_timer(15)=it_timer(15)+(j-i)

call system_clock(i,k)
CALL ENbond()
call system_clock(j,k)
it_timer(7)=it_timer(7)+(j-i)

!$omp section

call system_clock(i,k)
call NEIGHBORLIST(NMINCELL, NBUFFER, atype, pos)
call system_clock(j,k)
it_timer(5)=it_timer(5)+(j-i)

call system_clock(i,k)
CALL BOCALC(NMINCELL, NBUFFER, atype, pos)
call system_clock(j,k)
it_timer(6)=it_timer(6)+(j-i)

call system_clock(i,k)
CALL Ebond()
call system_clock(j,k)
it_timer(8)=it_timer(8)+(j-i)

call system_clock(i,k)
CALL Elnpr()
call system_clock(j,k)
it_timer(9)=it_timer(9)+(j-i)

call system_clock(i,k)
CALL E3b()
call system_clock(j,k)
it_timer(11)=it_timer(11)+(j-i)

call system_clock(i,k)
CALL E4b()
call system_clock(j,k)
it_timer(12)=it_timer(12)+(j-i)

call system_clock(i,k)
CALL Ehb()
call system_clock(j,k)
it_timer(10)=it_timer(10)+(j-i)

call system_clock(i,k)
CALL ForceBondedTerms(NMINCELL)
call system_clock(j,k)
it_timer(13)=it_timer(13)+(j-i)

!$omp end sections
!$omp end parallel

! merge force values from ENbond
f(:,:)=f(:,:)+fnb(:,:)

call system_clock(i,k)
dr(1:3)=0.d0
CALL COPYATOMS(MODE_CPBK,dr, NBUFFER, atype, pos, vdummy, f, q) 
call system_clock(j,k)
it_timer(14)=it_timer(14)+(j-i)

!rr(1:3)=0.d0
!do i=1, NATOMS
!   rr(1:3)=rr(1:3)+f(1:3,i)
!   print'(i,6f)',i,f(1:3,i),q(i)
!enddo
!print'(a,3f)','frcsum: ', rr(1:3)
!stop

!--- calculate kinetic part of stress components and add to <astr>.
#ifdef STRESS
call stress()
#endif

return

CONTAINS 

!----------------------------------------------------------------------
subroutine ForceBondedTerms(nlayer)
use atoms; use parameters
!----------------------------------------------------------------------
implicit none
integer,intent(IN) :: nlayer
integer :: c1,c2,c3 , i,i1,j,j1
real(8) :: dr(3), dr2, ff(3)

DO c1=-nlayer, cc(1)-1+nlayer
DO c2=-nlayer, cc(2)-1+nlayer
DO c3=-nlayer, cc(3)-1+nlayer

  i = header(c1, c2, c3)
  do i1=1, nacell(c1, c2, c3)

     do j1=1, nbrlist(i,0)
        j=nbrlist(i,j1)
        dr(1:3) = pos(1:3,i) - pos(1:3,j)
        ff(1:3) = ccbnd(i)*dBOp(i,j1)*dr(1:3)
        f(1:3,i) = f(1:3,i) - ff(1:3)
        f(1:3,j) = f(1:3,j) + ff(1:3)
#ifdef STRESS
        ia=i; ja=j
        include 'stress'
#endif
     enddo

!--- reset ccbnd to zero for next turn. first rest is done during
!initialization.
     ccbnd(i)=0.d0

     i = llist(i)
  enddo
enddo; enddo; enddo

end subroutine


!-------------------------------------------------------------------------------------
subroutine Elnpr()
use parameters;use atoms
!-------------------------------------------------------------------------------------
implicit none

integer :: c1,c2,c3,ii
integer :: i,i1,j,j1,j2,ity, jty, inxn, idEh
real(8) :: coeff(3)

!--- Lone Pair Energy Terms
real(8) :: Clp, CElp(1), PElp, dEh, mdEh
real(8) :: explp1,cnlp1,expvd2, dElp, deltaE
real(8),allocatable :: deltalp(:)

!--- Overcoordination Energy Terms
real(8) :: sum_ovun1, sum_ovun2
real(8) :: deltalpcorr, PEover, DlpV_i
real(8) :: expovun2, expovun1

!--- Undercoordination Energy Terms
real(8) :: expovun2n, expovun6, expovun8
real(8) :: PEunder

real(8) :: CEover(7), CEunder(6),fsum(1:3), CElp_b, CElp_d, CElp_bpp
integer :: n,n1,nty, m,m1,mty

real(8) :: div_expovun2, div_expovun2n, div_expovun1, div_expovun8

real(8) :: get_exp

allocate(deltalp(NBUFFER),stat=ast)

!=== preparation ==============================================================
do i = 1, NBUFFER
   ity = atype(i)

   if(ity==0) cycle

   deltaE = -Vale(ity) + Val(ity) + delta(i)

   dEh = deltaE*0.5d0
   idEh = is_idEH*int(dEh)
   explp1 = exp( -plp1(ity)*(2.d0 + deltaE - 2*idEh)**2 )
   Clp = 2.d0*plp1(ity)*explp1*(2.d0 + deltaE - 2*idEh)

   dDlp(i) = Clp

!--- delta function ver.
!     dDlp(i) = Clp + (0.5d0 - Clp)*( epsrn*abs(mdEh)**(epsrn-1.d0) )

   nlp(i) = explp1 - dble(idEh)
   deltalp(i) = nlpopt(ity) - nlp(i)

!--- if mass is greater than 21.0, deltalp(i) becomes zero. see poten.f line 750.
!--- if (amas(ity1).gt.21.0) dfvl=0.0d0
!--- diffvlph=dfvl*(voptlp-vlptemp(i1))

   if(mass(ity)>21.d0) deltalp(i) = 0.d0

enddo
!============================================================== preparation ===

PE(2:4)=0.d0
do i=1, NATOMS
   ity = atype(i)

   sum_ovun1 = 0.d0
   sum_ovun2 = 0.d0
   do j1 = 1, nbrlist(i,0)
      j = nbrlist(i, j1)
      jty = atype(j)
      inxn = inxn2(ity,jty)
      sum_ovun1 = sum_ovun1 + povun1(inxn)*Desig(inxn)*BO(0,i,j1)
      sum_ovun2 = sum_ovun2 + (delta(j)-deltalp(j))*(BO(2,i,j1) + BO(3,i,j1))
   enddo

!--- Lone Pair
   expvd2 = exp(-75.d0*deltalp(i))
   dElp = plp2(ity)*( (1.d0 + expvd2) + 75.d0*deltalp(i)*expvd2 ) / (1.d0 + expvd2)**2

!--- Over Coordinate + Common part with Under Coordinate
   expovun1 = povun3(ity) * exp( povun4(ity) * sum_ovun2 )
   deltalpcorr = delta(i) - deltalp(i)/(1.d0 + expovun1) 
   expovun2=exp(povun2(ity)*deltalpcorr)

!=== if one atom flys away, this term becomes zero because the total bond-order becomes zero.
!=== Add a small value in the denominator to avoid it. See poten.f line 787,
!=== hulpp=(1.0/(vov1+aval(ity1)+1e-8))
   DlpV_i = 1.d0/(deltalpcorr + Val(ity) + 1.d-8)

!--- Under Coordinate
   expovun2n = 1.d0 / expovun2
   expovun6 = exp(povun6(ity)*deltalpcorr)
   expovun8 = povun7(ity)*exp( povun8(ity)*sum_ovun2 )

   div_expovun1  = 1.d0/(1.d0 + expovun1)
   div_expovun2  = 1.d0/(1.d0 + expovun2)
   div_expovun2n = 1.d0/(1.d0 + expovun2n)
   div_expovun8  = 1.d0/(1.d0 + expovun8)

!--- Energy Calculation
   PElp = plp2(ity)*deltalp(i)/(1.d0+expvd2)
   PEover = sum_ovun1 * DlpV_i * deltalpcorr * div_expovun2
   PEunder = -povun5(ity) * (1.d0 - expovun6)*div_expovun2n*div_expovun8

!--- if the representitive atom is a resident, sum thier potential energies.
   PE(2) = PE(2) + PElp
   PE(3) = PE(3) + PEover
   PE(4) = PE(4) + PEunder

!--- Coefficient Calculation
   CElp(1) = dElp*dDlp(i)

   CEover(1) = deltalpcorr*DlpV_i*div_expovun2
   CEover(2) = sum_ovun1*DlpV_i*div_expovun2 * &
         (1.d0 - deltalpcorr*DlpV_i - povun2(ity)*deltalpcorr*div_expovun2n )

   CEover(3) = CEover(2)*( 1.d0 - dDlp(i)*div_expovun1 )
   CEover(4) = CEover(2)*deltalp(i)*povun4(ity)*expovun1*div_expovun1**2 

   CEunder(1) = (  povun5(ity)*povun6(ity)*expovun6*div_expovun8 &
                 + PEunder*povun2(ity)*expovun2n)*div_expovun2n 
   CEunder(2) =-PEunder*povun8(ity)*expovun8*div_expovun8
   CEunder(3) = CEunder(1)*(1.d0 - dDlp(i)*div_expovun1 )
   CEunder(4) = CEunder(1)*deltalp(i)*povun4(ity)*expovun1*div_expovun1**2 + CEunder(2)

!          CElp(:)=0.d0;    PE(2)=0.d0
!          CEover(:)=0.d0;  PE(3)=0.d0
!          CEunder(:)=0.d0; PE(4)=0.d0

!--- Force Calculation
   do j1 = 1, nbrlist(i,0) 
      j = nbrlist(i, j1) 
      jty = atype(j)
      inxn = inxn2(ity,jty)

      CEover(5) = CEover(1)*povun1(inxn)*desig(inxn)
      CEover(6) = CEover(4)*( 1.d0 - dDlp(j)) * (BO(2,i,j1) + BO(3,i,j1))
      CEover(7) = CEover(4)*( delta(j) - deltalp(j) )

      CEunder(5) = CEunder(4)*(1.d0 - dDlp(j))*(BO(2,i,j1) + BO(3,i,j1))
      CEunder(6) = CEunder(4)*(delta(j)-deltalp(j)) 
 
      CElp_b = CElp(1) + CEover(3) + CEover(5) + CEunder(3)
      CElp_bpp = CEover(7) + CEunder(6)
      coeff(1:3) = CElp_b + (/0.d0, CElp_bpp, CElp_bpp/) 

      i1=nbrindx(i,j1)
      call ForceBbo(i,j1, j,i1, coeff)

      CElp_d  = CEover(6) + CEunder(5)
      call ForceD(j, CElp_d)
   enddo

enddo ! i-loop

deallocate(deltalp,stat=ast)

END subroutine

!------------------------------------------------------------------------------------
subroutine E3b()
use parameters;use atoms
!------------------------------------------------------------------------------------
implicit none
integer :: i,j,k, i1,j1,k1, i2,j2,k2, ity,jty,kty, inxn, inxnhb, n,n1
integer :: ii,c1,c2,c3,NB

!--- Valency Energy Calculation Terms:
real(8) :: PEval, fn7ij, fn7jk, fn8j, delta_ang
real(8) :: rij(0:3), rjk(0:3), rij_mag, rjk_mag, SBO2, SBO 
real(8) :: sum_BO8, theta_ijk, theta0, theta_diff, exp2 
real(8) :: sin_ijk, cos_ijk 
real(8) :: Cf7ij, Cf7jk, Cf8j, Ctheta_diff, Ctheta0, CSBO2, dSBO(2)
real(8) :: BOij_p4, BOjk_p4, exp3ij, exp3jk, exp6, exp7, trm8
real(8) :: sum_SBO1, prod_SBO

!--- Penalty Energy Calculation Terms:
real(8) :: fn9, PEpen
real(8) :: exp_pen3, exp_pen4, exp_pen2ij, exp_pen2jk
real(8) :: trm_pen34, Cf9j

!--- Conjugation Energy (3body) Calculation Terms:
real(8) :: delta_val, PEcoa
real(8) :: sum_BOi, sum_BOk, exp_coa3i, exp_coa3k, exp_coa4i, exp_coa4k
real(8) :: exp_coa2 

!--- misc
real(8) :: CEval(9), CEpen(3), CEcoa(5)
real(8) :: explp1, Cnlp1, deltaE 
real(8) :: CE3body_d(3), CE3body_b(3), CE3body_bpp, CE3body_a, fsum(3), coeff(3)
real(8) :: BOij,BOjk

real(8) :: get_exp
integer :: l2g

!NOTICE: <cutof2> is used to get exactly the same energy in original reaxFF code. 
real(8) :: cutof2
cutof2 = cutof2_esub

PE(5:7)=0.d0
do j=1, NATOMS
   jty = atype(j)

   sum_BO8 = 0.d0
   sum_SBO1 = 0.d0
   do n1=1, nbrlist(j,0)
      sum_BO8 = sum_BO8 - BO(0,j,n1)**8.d0
      sum_SBO1 = sum_SBO1 + BO(2,j,n1) + BO(3,j,n1)
   enddo
   prod_SBO = exp(sum_BO8)

   delta_ang = delta(j) + Val(jty) - Valangle(jty)   

   do i1=1, nbrlist(j,0)-1

!--- NOTICE: Subtract the bond-order between i-j by cutof2 and use it for energy and force calc.
      BOij = BO(0,j,i1) - cutof2
      if(BOij>0.d0) then ! react.f, line 4827 
      i=nbrlist(j,i1)
      ity = atype(i)

      rij(1:3) = pos(1:3,i) - pos(1:3,j)
      rij(0) = sqrt( sum(rij(1:3)*rij(1:3)) )

      do k1=i1+1, nbrlist(j,0)

!--- NOTICE: Subtract the bond-order between i-j by cutof2 and use it for energy and force calc.
         BOjk = BO(0,j,k1)-cutof2

         if(BOjk>0.d0) then !react.f, line 4830
         if(BO(0,j,i1)*BO(0,j,k1)>cutof2) then !react.f, line 4831

         k=nbrlist(j,k1)
         kty = atype(k)

         rjk(1:3) = pos(1:3,j) - pos(1:3,k)
         rjk(0) = sqrt( sum(rjk(1:3)*rjk(1:3)) )  

         cos_ijk = -sum( rij(1:3)*rjk(1:3) ) / ( rij(0) * rjk(0) ) 
         if(cos_ijk>MAXANGLE) cos_ijk = MAXANGLE
         if(cos_ijk<MINANGLE) cos_ijk = MINANGLE
         theta_ijk = acos(cos_ijk) 
         sin_ijk = sin(theta_ijk)

!--- Check the type of 3atoms combination.
         inxn = inxn3(ity,jty,kty)
         if(inxn /= 0) then

!--- PEval part:
            BOij_p4 = BOij**pval4(inxn)
            exp3ij = exp( -pval3(jty)*BOij_p4 )
            fn7ij = 1.d0 - exp3ij
            BOjk_p4 = BOjk**pval4(inxn)
            exp3jk = exp( -pval3(jty)*BOjk_p4 )
            fn7jk = 1.d0 - exp3jk

            exp6 = exp( pval6(inxn)*delta_ang )
            exp7 = exp(-pval7(inxn)*delta_ang )
            trm8 = 1.d0 + exp6 + exp7
            fn8j = pval5(jty) - (pval5(jty)-1.d0)*(2.d0 + exp6)/trm8 

            SBO = sum_SBO1 + (1.d0 - prod_SBO)*(-delta_ang - pval8(inxn)*nlp(j))
            if(SBO.LE.0) SBO2 = 0.d0
            if(SBO.GT.0) SBO2 = SBO**pval9(inxn)
            if(SBO.GT.1) SBO2 = 2.d0 - (2.d0 - SBO)**pval9(inxn)
            if(SBO.GT.2) SBO2 = 2.d0
 
            theta0 = pi - theta00(inxn) * (1.d0 - exp( -pval10(inxn)*(2.d0-SBO2) ))
            theta_diff = theta0 - theta_ijk
            exp2 = exp( -pval2(inxn) * theta_diff*theta_diff ) 
   
            PEval = fn7ij * fn7jk * fn8j * (pval1(inxn) - pval1(inxn)*exp2)   
   
!---  PEval derivative part:
            Cf7ij = pval3(jty)*pval4(inxn)*( BOij**(pval4(inxn)-1.d0) )*exp3ij 
            Cf7jk = pval3(jty)*pval4(inxn)*( BOjk**(pval4(inxn)-1.d0) )*exp3jk

            Cf8j = (1.d0 - pval5(jty))/(trm8*trm8) *                                &
                     (                                                              &
                          pval6(inxn)*exp6*trm8                                     &
                        -(2.d0 + exp6)*(pval6(inxn)*exp6 - pval7(inxn)*exp7)        &
                      )
    
            Ctheta_diff = 2.d0*pval2(inxn)*theta_diff*exp2/(1.d0 - exp2)
            Ctheta0 = pval10(inxn)*theta00(inxn)*exp( -pval10(inxn)*(2.d0 - SBO2) )  

            if((SBO.LE.0).OR.(SBO.GT.2))  CSBO2 = 0.d0
            if((SBO.GT.0).AND.(SBO.LE.1)) CSBO2 = pval9(inxn)*SBO**(pval9(inxn)-1.d0)
            if((SBO.GT.1).AND.(SBO.LE.2)) CSBO2 = pval9(inxn)*(2.d0-SBO)**(pval9(inxn)-1.d0)
  
            dSBO(1) =-8.d0*prod_SBO*( delta_ang + pval8(inxn)*nlp(j) )
            dSBO(2) = (prod_SBO - 1.d0)*(1.d0 - pval8(inxn)*dDlp(j) )
 
            CEval(1) = Cf7ij * fn7jk * fn8j * pval1(inxn) * (1.d0 - exp2)
            CEval(2) = fn7ij * Cf7jk * fn8j * pval1(inxn) * (1.d0 - exp2) 
            CEval(3) = fn7ij * fn7jk * Cf8j * pval1(inxn) * (1.d0 - exp2)

            CEval(4) = 2.0d0*pval1(inxn)*pval2(inxn)*fn7ij*fn7jk*fn8j*exp2*theta_diff
            CEval(5) = CEval(4)*Ctheta0*CSBO2
            CEval(6) = CEval(5)*dSBO(1)
            CEval(7) = CEval(5)*dSBO(2)
            CEval(8) = CEval(4)/sin_ijk

!--- PEpen part
            exp_pen3 = exp(-ppen3(inxn)*delta(j))
            exp_pen4 = exp( ppen4(inxn)*delta(j))
            fn9 = (2.d0 + exp_pen3) / (1.d0 + exp_pen3 + exp_pen4)
            exp_pen2ij = exp( -ppen2(inxn)*(BOij-2.d0)*(BOij-2.d0) )
            exp_pen2jk = exp( -ppen2(inxn)*(BOjk-2.d0)*(BOjk-2.d0) )

            PEpen = ppen1(inxn) * fn9 * exp_pen2ij * exp_pen2jk 

!--- PEpen derivative part:
            trm_pen34 = 1.d0 + exp_pen3 + exp_pen4
            Cf9j = ( -ppen3(inxn)*exp_pen3*trm_pen34                                   &
                     -(2.d0 + exp_pen3)*(-ppen3(inxn)*exp_pen3 + ppen4(inxn)*exp_pen4) &
                   ) / (trm_pen34*trm_pen34)
            CEpen(1) = Cf9j / fn9 
            CEPen(2) = -2.d0 * ppen2(inxn) * (BOij - 2.d0)
            CEPen(3) = -2.d0 * ppen2(inxn) * (BOjk - 2.d0)
            CEpen(1:3) = CEpen(1:3)*PEpen 

!print'(4d25.13,3i5)', PEpen, CEPen(1:3), &
!         mod(l2g(atype(i))-1,168)+1, mod(l2g(atype(j))-1,168)+1,mod(l2g(atype(k))-1,168)+1

!--- PEcoa part:
            sum_BOi = delta(i) + Val(ity)
            sum_BOk = delta(k) + Val(kty)
            delta_val = delta(j) + Val(jty) - Valval(jty) 

            exp_coa2 = exp( pcoa2(inxn)*delta_val )
            exp_coa3i = exp( -pcoa3(inxn)*(-BOij + sum_BOi)**2 )
            exp_coa3k = exp( -pcoa3(inxn)*(-BOjk + sum_BOk)**2 )
            exp_coa4i = exp( -pcoa4(inxn)*( BOij - 1.5d0)**2 )
            exp_coa4k = exp( -pcoa4(inxn)*( BOjk - 1.5d0)**2 )

            PEcoa = pcoa1(inxn)/(1.d0 + exp_coa2)*exp_coa3i*exp_coa3k*exp_coa4i*exp_coa4k

!--- PEcoa derivative part: 
            CEcoa(1) =-2.d0*pcoa4(inxn)*(BOij-1.5d0) !dBOij
            CEcoa(2) =-2.d0*pcoa4(inxn)*(BOjk-1.5d0) !dBOjk
            CEcoa(3) =-pcoa2(inxn) * exp_coa2 / (1.d0 + exp_coa2) !dDj
            CEcoa(4) =-2.d0*pcoa3(inxn)*( -BOij + sum_BOi )
            CEcoa(5) =-2.d0*pcoa3(inxn)*( -BOjk + sum_BOk ) 
            CEcoa(1:5)=CEcoa(1:5)*PEcoa

!--- if the j-atom is a resident, count the potential energies.

            PE(5) = PE(5) + PEval
            PE(6) = PE(6) + PEpen
            PE(7) = PE(7) + PEcoa

!                  CEval(:)=0.d0; PE(5)=0.d0
!                  CEpen(:)=0.d0; PE(6)=0.d0
!                  CEcoa(:)=0.d0; PE(7)=0.d0


            CE3body_b(1) = CEpen(2) + CEcoa(1) - CEcoa(4) + CEval(1) !BO_ij 
            CE3body_b(2) = CEpen(3) + CEcoa(2) - CEcoa(5) + CEval(2) !BO_jk 

            CE3body_d(1) = CEpen(1) + CEcoa(3) + CEval(3) + CEval(7) !delta_j 

            CE3body_d(2) = CEcoa(4) !delta_i
            CE3body_d(3) = CEcoa(5) !delta_k

            CE3body_a = CEval(8)

!--- Force calculation
            j1 = nbrindx(j, i1)
            call ForceB(i,j1, j,i1, CE3body_b(1)) !BO_ij

            j1 = nbrindx(j, k1)
            call ForceB(j,k1, k,j1, CE3body_b(2)) !BO_jk

            do n1=1, nbrlist(j,0)
               coeff(1:3) = CE3body_d(1) + CEval(6)*BO(0,j,n1)**7 + (/0.d0, CEval(5),  CEval(5)/)

               n  = nbrlist(j, n1)
               j1 = nbrindx(j, n1)
               call ForceBbo(j,n1, n,j1, coeff) 
            enddo

            call ForceD(i, CE3body_d(2))
            call ForceD(k, CE3body_d(3))

            call ForceA3(CE3body_a, i, j, k, rij, rjk) 

         endif !inxn3
         endif ! if(BOij*BOjk>MINBO0) then
         endif ! if(BOjk>MINBO0) then
      enddo ! k-loop
      endif ! if(BOij>MINBO0) then
   enddo ! i-loop
enddo ! j-loop


END subroutine
!----------------------------------------------------------------------------------------------------------------------
subroutine Ehb()
use parameters;use atoms
! Note: 02-09-05 <kn>
! To find out hydrogen bonding combinations, <vnhbp> (one atom parameter, 2nd row, 8th column) is used to 
! identify whether an atoms is hydrogen or not, <vnhbp>=1 for H, <vnhbp>=2 for O,N,S and <vnhbp>=0 for others. 
! If a 3atom combination doesn't consists of H,O,N,S, the hydrogen bonding doens't exist.
! Among the three atoms, the center atom is always a hydrogen. In the original code, atom IDs are given 
! in such as j2 <-> j1(H) <=> i2. <-> is bonding interaction, <=> is non-bonding interaction.
! If the magnitude of the donor/acceptor bond is less than 1.d-2, the hydrogen bonding doesn't exit.
! If the interatomic distance between j2 <-> i2 is larger than 10[A], the interaction doesn't exit. 
!----------------------------------------------------------------------------------------------------------------------
implicit none

! <rchb>    cutoff length for atom j2 and i2 which is 10[A]

integer :: i,j,k, i1,j1,k1, ii,jj,kk, c1,c2,c3,c4,c5,c6, mn
integer :: ity,jty,kty,inxnhb
real(8) :: rij(0:3), rjk(0:3), rik(3), rik2
real(8) :: theta_ijk, cos_ijk, sin_ijk_half
real(8) :: sin_ijk, cos_xhz1, sin_xhz4, exp_hb2, exp_hb3
real(8) :: PEhb, CEhb(3), ff(3), dr(3)
real(8) :: get_exp

integer :: l2g
logical :: isEhb

PE(10)=0.d0
do c1=0, cc(1)-1
do c2=0, cc(2)-1
do c3=0, cc(3)-1

   i=header(c1,c2,c3)
   do ii=1,nacell(c1,c2,c3)
      ity = atype(i)

      do j1=1, nbrlist(i,0)
         j = nbrlist(i,j1)
         jty = atype(j)

         if( (jty==2) .and. (BO(0,i,j1)>MINBO0) ) then

             do kk=1, nbplist(i,0)

               k = nbplist(i,kk)

               kty = atype(k)

               inxnhb = inxn3hb(ity, jty, kty)

               if ( (j/=k).and.(i/=k).and.(inxnhb/=0) ) then

                  rik(1:3) = pos(1:3,i) - pos(1:3,k)
                  rik2 = sum(rik(1:3)*rik(1:3))

                  if(rik2<rchb2) then
      
                     rjk(1:3) = pos(1:3,j) - pos(1:3,k)
                     rjk(0) = sqrt( sum(rjk(1:3)*rjk(1:3)) )
      
                     rij(1:3) = pos(1:3,i) - pos(1:3,j)
                     rij(0) = sqrt( sum(rij(1:3)*rij(1:3)) )  
      
                     cos_ijk = -sum( rij(1:3)*rjk(1:3) ) / (rij(0) * rjk(0) ) 
                     if(cos_ijk>MAXANGLE) cos_ijk = MAXANGLE
                     if(cos_ijk<MINANGLE) cos_ijk = MINANGLE

                     theta_ijk = acos(cos_ijk)
      
                     sin_ijk_half = sin(0.5d0*theta_ijk)
                     sin_xhz4 = sin_ijk_half**4
                     cos_xhz1 = ( 1.d0 - cos_ijk )
      
                     exp_hb2 = exp( -phb2(inxnhb)*BO(0,i,j1) )
                     exp_hb3 = exp( -phb3(inxnhb)*(r0hb(inxnhb)/rjk(0) + rjk(0)/r0hb(inxnhb) - 2.d0) )
      
                     PEhb = phb1(inxnhb)*(1.d0 - exp_hb2)*exp_hb3*sin_xhz4

                     PE(10) = PE(10) + PEhb
      
                     CEhb(1) = phb1(inxnhb)*phb2(inxnhb)*exp_hb2*exp_hb3*sin_xhz4
                     CEhb(2) =-0.5d0*phb1(inxnhb)*(1.d0 - exp_hb2)*exp_hb3*cos_xhz1
                     CEhb(3) =-PEhb*phb3(inxnhb)*( -r0hb(inxnhb)/rjk(0)**2 + 1.d0/r0hb(inxnhb) )*(1.d0/rjk(0))

                     i1 = nbrindx(i,j1)
                     call ForceB(i,j1, j,i1, CEhb(1))
                     call ForceA3(CEhb(2), i,j,k, rij, rjk)
      
                     ff(1:3) = CEhb(3)*rjk(1:3)
                     f(1:3,j) = f(1:3,j) - ff(1:3)
                     f(1:3,k) = f(1:3,k) + ff(1:3)

!--- stress calculation
#ifdef STRESS
                     ia=j; ja=k; dr(1:3)=rjk(1:3)
                     include 'stress'
#endif

                  endif ! if(rik2<rchb2)
               endif

           enddo !do kk=1, nbplist(i,0)

         endif ! if(BO(0,j,i1)>MINBO0)
      enddo 

      i=llist(i)
   enddo
enddo; enddo; enddo

end subroutine

!----------------------------------------------------------------------------------------------------------
subroutine ENbond()
use parameters;use atoms
!----------------------------------------------------------------------------------------------------------
!  This subroutine calculates the energy and the forces due to the Van der Waals and Coulomb terms 
!----------------------------------------------------------------------------------------------------------
implicit none
integer :: i, j,j1, ity,jty,nn 
integer :: c1,c2,c3,c4,c5,c6, m,n, mn

real(8) :: PEvdw, PEclmb, ff(3)
real(8) :: dr(0:3), dr2, fsum(3)
real(8) :: CEvdw, CEclmb,qij

integer :: ii,l2g, iid, jid

integer :: inxn, itb, itb1
real(8) :: drtb, drtb1

PE(11:13)=0.d0
do c1=0, cc(1)-1
do c2=0, cc(2)-1
do c3=0, cc(3)-1

   i = header(c1,c2,c3)
   do m = 1, nacell(c1,c2,c3)
      iid = l2g(atype(i))
      ity = atype(i) 
      
      PE(13) = PE(13) + CEchrge*(chi(ity)*q(i) + 0.5d0*eta(ity)*q(i)**2)

       do j1 = 1, nbplist(i,0) 
            j = nbplist(i,j1)

            jid = l2g(atype(j))

            if(jid<iid) then

               dr(1:3) = pos(1:3,i) - pos(1:3,j)
               dr2 = sum(dr(1:3)*dr(1:3))

               if(dr2<=rctap2) then

                  jty=atype(j)

                  inxn = inxn2(ity, jty)
!                  dr(0) = sqrt(dr2)

!--- get table index and residual value
!                  itb = int(dr(0)*UDRi)
                  itb = int(dr2*UDRi)
                  itb1 = itb+1
                  drtb = dr2 - itb*UDR
                  drtb = drtb*UDRi
                  drtb1= 1.d0-drtb

!--- van del Waals:
                  PEvdw  = drtb1*TBL_Evdw(0,itb,inxn)  + drtb*TBL_Evdw(0,itb1,inxn)
                  CEvdw  = drtb1*TBL_Evdw(1,itb,inxn)  + drtb*TBL_Evdw(1,itb1,inxn)
!--- Coulomb:
                  qij = q(i)*q(j)
                  PEclmb = drtb1*TBL_Eclmb(0,itb,inxn) + drtb*TBL_Eclmb(0,itb1,inxn)
                  PEclmb = PEclmb*qij
                  CEclmb = drtb1*TBL_Eclmb(1,itb,inxn) + drtb*TBL_Eclmb(1,itb1,inxn)
                  CEclmb = CEclmb*qij

                  PE(11) = PE(11) + PEvdw
                  PE(12) = PE(12) + PEclmb

                  ff(1:3) = (CEvdw+CEclmb)*dr(1:3)
       
                  fnb(1:3,i) = fnb(1:3,i) - ff(1:3)
                  fnb(1:3,j) = fnb(1:3,j) + ff(1:3)

!--- stress calculation
#ifdef STRESS
                   ia=i; ja=j
                   include 'stress'
#endif

               endif

            endif

       enddo  !do j1 = 1, nbplist(i,0) 

      i=llist(i)
   enddo
enddo; enddo; enddo

END subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine Ebond()
use atoms; use parameters 
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,j, i1,j1, ity, jty, inxn
integer :: c1,c2,c3, n
real(8) :: exp_be12,  CEbo, PEbo, coeff(3), ff(3)
real(8) :: get_exp

integer :: l2g,iid,jid

PE(1)=0.d0
do i=1, NATOMS
   ity = atype(i)

   iid=l2g(atype(i))
   do j1 = 1, nbrlist(i,0)

      j = nbrlist(i,j1)
      jid=l2g(atype(j))
      if(jid<iid) then

        jty = atype(j)
        inxn = inxn2(ity, jty)

        exp_be12 = exp( pbe1(inxn)*( 1.d0 - BO(1,i,j1)**pbe2(inxn) ) )

        PEbo = - Desig(inxn)*BO(1,i,j1)*exp_be12 - Depi(inxn)*BO(2,i,j1) - Depipi(inxn)*BO(3,i,j1) 

        PE(1) = PE(1) + PEbo  !<kn>

        CEbo = -Desig(inxn)*exp_be12*( 1.d0 - pbe1(inxn)*pbe2(inxn)*BO(1,i,j1)**pbe2(inxn) )
        coeff(1:3)= (/ CEbo, -Depi(inxn), -Depipi(inxn) /)

        i1 = nbrindx(i,j1)
        call ForceBbo(i,j1, j,i1, coeff)

      endif

   enddo
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine E4b()
use atoms; use parameters
!--------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,j,k,l, i1,j1,k1,l1, k2, ity,jty,kty,lty, inxn, m, m1

!--- angles
real(8) :: cos_ijkl(3), cos_ijkl_sqr, cos_2ijkl, sin_ijkl, sin_ijk, sin_jkl, tan_ijk_i, tan_jkl_i
real(8) :: cos_ijk, cos_jkl, theta_ijk, theta_jkl, omega_ijkl

!--- vectors
real(8) :: rij(0:3), rjk(0:3), rkl(0:3)
real(8) :: crs_ijk(0:3), crs_jkl(0:3)

real(8) :: delta_ang_jk, delta_ang_j, delta_ang_k
real(8) :: exp_tor1, exp_tor3, exp_tor4, exp_tor34_i, fn10, fn11, dfn11, fn12, PEtors, PEconj, cmn
real(8) :: exp_tor2(3)

!--- coefficients
real(8) :: CEtors(9), Cconj, CEconj(6), C4body_a(3), C4body_b(3), C4body_b_jk(3)
real(8) :: coeff(3)

real(8) :: btb2 !! 
real(8) :: BOij, BOjk, BOkl
real(8) :: cutof2

real(8) :: get_exp

integer :: jid,kid,l2g
real(8) :: esft_e4, e_e4, drtb_e4, drtb1_e4 
integer :: itb_e4,itb1_e4

integer :: c1,c2,c3,mn

cutof2 = cutof2_esub

PE(8:9)=0.d0
do j=1,NATOMS

  jty = atype(j)
  delta_ang_j = delta(j) + Val(jty) - Valangle(jty)
  jid=l2g(atype(j))

  do k1=1, nbrlist(j,0) 

!--- NOTICE: Subtract the bond-order between i-j by cutof2 and use it for energy and force calc.
     BOjk = BO(0,j,k1) - cutof2 !<kn>
     if(BO(0,j,k1) > cutof2) then         !poten.f,line 1829,1830

        k=nbrlist(j,k1)
        kid=l2g(atype(k))

        if (jid<kid) then

        kty = atype(k)

        delta_ang_k = delta(k) + Val(kty) - Valangle(kty)
        delta_ang_jk = delta_ang_j + delta_ang_k  

        rjk(1:3) = pos(1:3,j) - pos(1:3,k)
        rjk(0) = sqrt( sum(rjk(1:3)*rjk(1:3)) )
         
        do i1=1, nbrlist(j,0)

!--- NOTICE: Subtract the bond-order between i-j by cutof2 and use it for energy and force calc.
           BOij = BO(0,j,i1) - cutof2 !<kn>

!--- NOTICE: cutoff condition to ignore bonding.
           if((BO(0,j,i1)>cutof2) .and. ((BO(0,j,i1)*BO(0,j,k1))>cutof2)) then !poten.f from iv() calculataion

           i=nbrlist(j,i1)

           if (i/=k) then

              ity = atype(i)

              rij(1:3) = pos(1:3,i) - pos(1:3,j)
              rij(0) = sqrt( sum(rij(1:3)*rij(1:3)) )

!--- Calculate the angle i-j-k
              cos_ijk =-sum(rij(1:3)*rjk(1:3))/(rij(0)*rjk(0))
              if(cos_ijk>MAXANGLE) cos_ijk = MAXANGLE
              if(cos_ijk<MINANGLE) cos_ijk = MINANGLE

              theta_ijk = acos( cos_ijk ) 
              sin_ijk  = sin(theta_ijk)
              tan_ijk_i = 1.d0/tan(theta_ijk)

              call cross_product(rij, rjk, crs_ijk)

              do l1=1, nbrlist(k,0)

!--- NOTICE: Subtract the bond-order between k-l by cutof2 and use it for energy and force calc.
                 BOkl = BO(0,k,l1) - cutof2 !<kn>

!--- NOTICE: cutoff condition to ignore bonding.
                 if((BO(0,k,l1)>cutof2).and.(BO(0,j,k1)*BO(0,k,l1)>cutof2)) then !poten.f,line 1829,1830

                 l=nbrlist(k,l1)
                 lty = atype(l)
                 inxn = inxn4(ity,jty,kty,lty)

                 if ((inxn/=0).and.(i/=l).and.(j/=l)) then

!--- NOTICE: cutoff condition to ignore bonding.
                 if( (BO(0,j,i1)*(BO(0,j,k1)**2)*BO(0,k,l1)) > MINBO0) then

                 rkl(1:3) = pos(1:3,k) - pos(1:3,l)
                 rkl(0) = sqrt( sum(rkl(1:3)*rkl(1:3)) )

                 exp_tor2(1) = exp(-ptor2(inxn)*BOij)  ! i-j
                 exp_tor2(2) = exp(-ptor2(inxn)*BOjk)  ! j-k
                 exp_tor2(3) = exp(-ptor2(inxn)*BOkl)  ! k-l

                 exp_tor3 = exp(-ptor3(inxn)*delta_ang_jk )
                 exp_tor4 = exp( ptor4(inxn)*delta_ang_jk )
                 exp_tor34_i = 1.d0/(1.d0 + exp_tor3 + exp_tor4)

                 fn10 = (1.d0 - exp_tor2(1))*(1.d0 - exp_tor2(2))*(1.d0 - exp_tor2(3))
                 fn11 = (2.d0 + exp_tor3)/(1.d0 + exp_tor3 + exp_tor4)

                 fn12 = exp(-pcot2(inxn)*( (BOij-1.5d0)**2 + &
                                           (BOjk-1.5d0)**2 + &
                                           (BOkl-1.5d0)**2 ) )

!--- NOTICE: pi-bond value used here is not the subtracted one but the original value. 
                 btb2 = 2.d0 - BO(2,j,k1) - fn11         !<kn>
                 exp_tor1 = exp( ptor1(inxn)*btb2**2 )   !<kn>


!--- Get angle variables i-j-k, j-k-l, i-j-k-l
                 cos_jkl =-sum(rjk(1:3)*rkl(1:3))/(rjk(0)*rkl(0))
                 if(cos_jkl>MAXANGLE) cos_jkl = MAXANGLE
                 if(cos_jkl<MINANGLE) cos_jkl = MINANGLE

                 theta_jkl = acos( cos_jkl )
                 sin_jkl = sin(theta_jkl)
                 tan_jkl_i = 1.d0/tan(theta_jkl)

                 call cross_product(rjk, rkl, crs_jkl)

                 cos_ijkl(1) = sum( crs_ijk(1:3)*crs_jkl(1:3) )/( crs_ijk(0)*crs_jkl(0) )
                 if(cos_ijkl(1)>MAXANGLE) cos_ijkl(1) = MAXANGLE
                 if(cos_ijkl(1)<MINANGLE) cos_ijkl(1) = MINANGLE

                 omega_ijkl = acos(cos_ijkl(1))
                 cos_ijkl_sqr = cos_ijkl(1)*cos_ijkl(1)
                 cos_2ijkl = cos(2.d0*omega_ijkl)
                 cos_ijkl(2) = 1.d0 - cos_2ijkl
                 cos_ijkl(3) = 1.d0 + cos(3.d0*omega_ijkl)
                 sin_ijkl = sin(omega_ijkl)

                 PEtors = 0.5d0*fn10*sin_ijk*sin_jkl* &
                   ( V1(inxn)*(1.d0 + cos_ijkl(1) ) + V2(inxn)*exp_tor1*cos_ijkl(2) + V3(inxn)*cos_ijkl(3) )  !<kn>

                 PEconj = pcot1(inxn)*fn12*(1.d0 + (cos_ijkl_sqr - 1.d0)*sin_ijk*sin_jkl)

                 PE(8) = PE(8) + PEtors
                 PE(9) = PE(9) + PEconj

!--- Force coefficient calculation
!--- Torsional term
                 CEtors(1) = 0.5d0*sin_ijk*sin_jkl*(                                  &
                                                      V1(inxn)*(1.d0 + cos_ijkl(1))   & 
                                                    + V2(inxn)*exp_tor1*cos_ijkl(2)   &
                                                    + V3(inxn)*cos_ijkl(3)            &
                                                   )

                 CEtors(2) =-ptor1(inxn)*fn10*sin_ijk*sin_jkl*V2(inxn)*exp_tor1*btb2*cos_ijkl(2)

                 dfn11 = (-ptor3(inxn)*exp_tor3 + &
                     (ptor3(inxn)*exp_tor3 - ptor4(inxn)*exp_tor4)*(2.d0 + exp_tor3)*exp_tor34_i )*exp_tor34_i

                 CEtors(3) = CEtors(2)*dfn11

                 CEtors(4) = CEtors(1)*ptor2(inxn)*exp_tor2(1)*(1.d0 - exp_tor2(2))*(1.d0 - exp_tor2(3))
                 CEtors(5) = CEtors(1)*ptor2(inxn)*(1.d0 - exp_tor2(1))*exp_tor2(2)*(1.d0 - exp_tor2(3))
                 CEtors(6) = CEtors(1)*ptor2(inxn)*(1.d0 - exp_tor2(1))*(1.d0 - exp_tor2(2))*exp_tor2(3)

                 cmn = -0.5d0*fn10*( V1(inxn)*(1.d0 + cos_ijkl(1)) + V2(inxn)*exp_tor1*cos_ijkl(2) + V3(inxn)*cos_ijkl(3) )

                 CEtors(7) = cmn*sin_jkl*tan_ijk_i
                 CEtors(8) = cmn*sin_ijk*tan_jkl_i
                 CEtors(9) = fn10*sin_ijk*sin_jkl * &
                   ( 0.5d0*V1(inxn) - 2.d0*V2(inxn)*exp_tor1*cos_ijkl(1) + 1.5d0*V3(inxn)*(cos_2ijkl + 2.d0*cos_ijkl_sqr) )

!--- Conjugation Energy
                 Cconj = -2.d0*pcot2(inxn)*PEconj
                 CEconj(1) = Cconj*( BOij-1.5d0 )
                 CEconj(2) = Cconj*( BOjk-1.5d0 )
                 CEconj(3) = Cconj*( BOkl-1.5d0 )

                 CEconj(4) =-pcot1(inxn)*fn12*(cos_ijkl_sqr-1.d0)*tan_ijk_i*sin_jkl
                 CEconj(5) =-pcot1(inxn)*fn12*(cos_ijkl_sqr-1.d0)*sin_ijk*tan_jkl_i
                 CEconj(6) = 2.d0*pcot1(inxn)*fn12*cos_ijkl(1)*sin_ijk*sin_jkl

!              CEtors(:)=0.d0; PE(8)=0.d0
!              CEconj(:)=0.d0; PE(9)=0.d0

                 C4body_b(1:3) = CEconj(1:3) + CEtors(4:6) !dBOij, dBOjk, dBOkl
                 C4body_a(1:3) = CEconj(4:6) + CEtors(7:9) !ijk, jkl, ijkl

                 call ForceD(j, CEtors(3))
                 call ForceD(k, CEtors(3))

                 j1 = nbrindx(j, i1)
                 call ForceB(i,j1, j,i1, C4body_b(1))

!--- To take care of the derivative of BOpi(j,k), add <Ctors(2)> to 
!--- the full BOjk derivative coefficient <C4body_b(2)>, but only pi-bond 
!--- component.
                 C4body_b_jk(1:3) = C4body_b(2) + (/0.d0, CEtors(2), 0.d0 /)

                 j1 = nbrindx(j, k1)
                 call ForceBbo(j,k1, k,j1, C4body_b_jk)

                 k2 = nbrindx(k, l1)
                 call ForceB(k,l1, l,k2, C4body_b(3))

                 call ForceA3(C4body_a(1), i,j,k, rij,rjk)
                 call ForceA3(C4body_a(2), j,k,l, rjk,rkl)
                 call ForceA4(C4body_a(3), i,j,k,l, rij,rjk,rkl)

!              if(icheck(0)) call checker(i, j, k, l, PEtors, f(1:3,j) )

                 endif ! if( (BO(0,j,i1)*BO(0,j,k1)**2*BO(0,k,l1)) >MINBO0)
                 endif ! if ((inxn/=0).and.(i/=l).and.(j/=l))
                 endif ! if(BOkl>MINBO0)

               enddo ! l-loop

               endif ! if (i/=k) 
            endif ! if(BOij>MINBO0)

         enddo !i-loop

         endif !if(BOjk>MINBO0)
      endif !if(j>k)

   enddo ! k-loop

enddo

end subroutine

!-----------------------------------------------------------------------------------------
subroutine ForceD(i, coeff)
use atoms
! Calculate force from derivative of delta(i). 
!-----------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: i
real(8),intent(IN) :: coeff
integer :: i1, j,j1, n,n1, inbrs, jnbrs
real(8) :: Cbond(3), dr(3), ff(3)

do j1=1, nbrlist(i,0)
  j  = nbrlist(i,j1)
  i1 = nbrindx(i,j1)

  Cbond(1) = coeff*(A0(i,j1) + BO(0,i,j1)*A1(i,j1) )! Coeff of BOp
  dr(1:3) = pos(1:3,i)-pos(1:3,j)
  ff(1:3) = Cbond(1)*dBOp(i,j1)*dr(1:3)
  f(1:3,i) = f(1:3,i) - ff(1:3)
  f(1:3,j) = f(1:3,j) + ff(1:3)

#ifdef STRESS
  ia=i; ja=j
  include 'stress'
#endif

  Cbond(2)=coeff*BO(0,i,j1)*A2(i,j1) ! Coeff of deltap_i
  Cbond(3)=coeff*BO(0,i,j1)*A2(j,i1) ! Coeff of deltap_j

  ccbnd(i)=ccbnd(i)+Cbond(2)
  ccbnd(j)=ccbnd(j)+Cbond(3)
!--- new ---

enddo

return
end subroutine

!-----------------------------------------------------------------------------------------
subroutine ForceB(i,j1, j,i1,coeff)
use atoms
! Derivative of BOij using new bond order definition. Only difference is that sigma BO
! prime is replaced with full BOp. The derivative of BO becomes a bit simpler due to 
! the new definition.
!-----------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: i,j1, j,i1
real(8),intent(IN) :: coeff
integer :: n, n1
real(8) :: Cbond(3),dr(3),ff(3)
! optimization
real(8) :: ri(3),rj(3),fi(3),fj(3)

Cbond(1) = coeff*(A0(i,j1) + BO(0,i,j1)*A1(i,j1) )! Coeff of BOp

dr(1:3) = pos(1:3,i) - pos(1:3,j)
ff(1:3) = Cbond(1)*dBOp(i,j1)*dr(1:3)
f(1:3,i) = f(1:3,i) - ff(1:3)
f(1:3,j) = f(1:3,j) + ff(1:3)

#ifdef STRESS
ia=i; ja=j
include 'stress'
#endif

!--- A3 is not necessary anymore with the new BO def. 
Cbond(2)=coeff*BO(0,i,j1)*A2(i,j1) ! Coeff of deltap_i
Cbond(3)=coeff*BO(0,i,j1)*A2(j,i1) ! Coeff of deltap_j

ccbnd(i)=ccbnd(i)+Cbond(2)
ccbnd(j)=ccbnd(j)+Cbond(3)

return
end subroutine

!-----------------------------------------------------------------------------------------
subroutine ForceBbo(i,j1, j,i1, coeff)
use atoms; use parameters
! Calculate force from derivative of BOij using different coefficient values 
!-----------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: i,j1, j,i1
real(8),intent(IN) :: coeff(3)
integer :: n,n1
real(8) :: Cbond(3),dr(3), ff(3),cBO(3),cf(3)

! optimization
real(8) :: ri(3),rj(3),fi(3),fj(3)

!--- With the new bond-order definition, 1st term is the derivative of
!"full"-bond order,
!--- 2nd is for pi-bond order and 3rd is for pipi-bond order.
cf(1:3) = (/ coeff(1), coeff(2)-coeff(1), coeff(3)-coeff(1) /)

Cbond(1) = cf(1)*(A0(i,j1) + BO(0,i,j1)*A1(i,j1))*dBOp(i,j1)        & !full BO
         + cf(2)*BO(2,i,j1)*( dln_BOp(2,i,j1)+A1(i,j1)*dBOp(i,j1) ) & !pi   BO
         + cf(3)*BO(3,i,j1)*( dln_BOp(3,i,j1)+A1(i,j1)*dBOp(i,j1) )   !pipi BO

dr(1:3) = pos(1:3,i)-pos(1:3,j)
ff(1:3) = Cbond(1)*dr(1:3)

f(1:3,i) = f(1:3,i) - ff(1:3)
f(1:3,j) = f(1:3,j) + ff(1:3)

#ifdef STRESS
ia=i; ja=j
include 'stress'
#endif

!--- 1st element is "full"-bond order.
cBO(1:3) = (/cf(1)*BO(0,i,j1),  cf(2)*BO(2,i,j1),  cf(3)*BO(3,i,j1) /)

Cbond(2)=cBO(1)*A2(i,j1) + (cBO(2)+cBO(3))*A3(i,j1)
Cbond(3)=cBO(1)*A2(j,i1) + (cBO(2)+cBO(3))*A3(j,i1)

ccbnd(i)=ccbnd(i)+Cbond(2)
ccbnd(j)=ccbnd(j)+Cbond(3)

return
end subroutine


!-----------------------------------------------------------------------------------------
subroutine ForceA4(coeff, i, j, k, l, da0, da1, da2)
use atoms
! derivative of <cos_ijkl>
!-----------------------------------------------------------------------------------------
implicit none
! Daa(0)=Daa-1, Daa(-1)=Da-1a-2
! Caa( 0,0)=Caa,   Caa( 0,-1)=Caa-1,   Caa( 0,-2)=Caa-2,
! Caa(-1,0)=Ca-1a, Caa(-1,-1)=Ca-1a-1, Caa(-1,-2)=Ca-1a-2
! Caa(-2,0)=Ca-2a, Caa(-2,-1)=Ca-2a-1, Caa(-2,-2)=Ca-2a-2
! da0 = ri - rj,  da1 = rj - rk,  da2 = rk - rl

integer,INTENT(IN) :: i,j,k,l
real(8),INTENT(IN) :: coeff, da0(0:3), da1(0:3), da2(0:3)
real(8) :: Daa(-1:0), Caa(-2:0,-2:0),Cwi(3), Cwj(3), Cwk(3), Cwl(3)
real(8) :: DDisqr, coDD, com
real(8) :: fij(3), fik(3), fil(3), fjk(3), fjl(3), fkl(3)
real(8) :: rij(3), rik(3), ril(3), rjk(3), rjl(3), rkl(3)
real(8) :: dr(3), ff(3)

Caa( 0,0)=da0(0)*da0(0);Caa( 0,-1)=sum(da0(1:3)*da1(1:3));Caa( 0,-2)=sum(da0(1:3)*da2(1:3))
Caa(-1,0)=Caa(0,-1);    Caa(-1,-1)=da1(0)*da1(0);         Caa(-1,-2)=sum(da1(1:3)*da2(1:3))
Caa(-2,0)=Caa(0,-2);    Caa(-2,-1)=Caa(-1,-2);            Caa(-2,-2)=da2(0)*da2(0)

Daa( 0) = Caa( 0, 0)*Caa(-1,-1) - Caa( 0,-1)*Caa( 0,-1)
Daa(-1) = Caa(-1,-1)*Caa(-2,-2) - Caa(-1,-2)*Caa(-1,-2)

DDisqr = 1.d0/sqrt( Daa(0)*Daa(-1) )
coDD = coeff*DDisqr

com = Caa(-1, 0)*Caa(-1,-2) - Caa(0,-2)*Caa(-1,-1)

!--- Some of calculations are unnecessary due to the action-reaction relation.
Cwi(1) = Caa(-1,-1)  / Daa(0)*com
Cwi(2) =-(Caa(-1,-2) + Caa(0,-1)/Daa(0)*com)
Cwi(3) = Caa(-1,-1)

Cwj(1) =-( Caa(-1,-2) + (Caa(-1,-1) + Caa(-1,0))/Daa(0)*com )
Cwj(2) =-(-Caa(-1,-2) - 2*Caa(0,-2) - Caa(-2,-2)/Daa(-1)*com - (Caa(0,0)+Caa(-1,0))/Daa(0)*com )
Cwj(3) =-( Caa(-1, 0) + Caa(-1,-1)  + Caa(-1,-2)/Daa(-1)*com )

!Cwk(1) = ( Caa(-2,-1) + Caa(-1,-1)  + Caa(-1,0)/Daa(0)*com )
!Cwk(2) =-( Caa(-1, 0) + 2*Caa(-2,0) + (Caa(-2,-2) + Caa(-2,-1))/Daa(-1)*com + Caa(0,0)/Daa(0)*com )
!Cwk(3) = ( Caa(-1, 0) + (Caa(-1,-1) + Caa(-2,-1))/Daa(-1)*com )

Cwl(1) =-Caa(-1,-1)
Cwl(2) = ( Caa(-1,0)  + Caa(-2,-1)/Daa(-1)*com )
Cwl(3) =-( Caa(-1,-1) / Daa(-1)*com )

rij(1:3) = da0(1:3)
rjk(1:3) = da1(1:3)
rkl(1:3) = da2(1:3)

fij(1:3) = coDD*( Cwi(1)*rij(1:3) + Cwi(2)*rjk(1:3) + Cwi(3)*rkl(1:3) )
fjk(1:3) = coDD*( (Cwj(1)+Cwi(1))*rij(1:3) + (Cwj(2)+Cwi(2))*rjk(1:3) + (Cwj(3)+Cwi(3))*rkl(1:3) )
fkl(1:3) =-coDD*( Cwl(1)*rij(1:3) + Cwl(2)*rjk(1:3) + Cwl(3)*rkl(1:3) ) 

f(1:3,i) = f(1:3,i) + fij(1:3) 
f(1:3,j) = f(1:3,j) - fij(1:3) + fjk(1:3) 
f(1:3,k) = f(1:3,k) - fjk(1:3) + fkl(1:3)
f(1:3,l) = f(1:3,l) - fkl(1:3)

!--- stress calculation
#ifdef STRESS
ia=i; ja=j; dr(1:3)=rij(1:3); ff(1:3)=-fij(1:3)
include 'stress'
ia=j; ja=k; dr(1:3)=rjk(1:3); ff(1:3)=-fjk(1:3)
include 'stress'
ia=k; ja=l; dr(1:3)=rkl(1:3); ff(1:3)=-fkl(1:3)
include 'stress'
#endif

!--- Check N3rd ---
!  print'(a,5f20.13)','N3rd: ',Cwi(1)-Cwi(2)+Cwj(1), Cwi(2)-Cwi(3)+Cwk(1), &
!        Cwi(3)+Cwl(1), Cwj(2)-Cwj(3)-Cwk(1)+Cwk(2), Cwk(3)-Cwl(2)+Cwl(3)

end subroutine

!-----------------------------------------------------------------------
subroutine ForceA3(coeff,i,j,k,da0, da1)
use atoms
! derivative of <cos_ijk>
!-----------------------------------------------------------------------
implicit none
! Caa(0,0)=Caa, Caa(0,-1)=Caa-1, Caa(-1,0)=Ca-1a, Caa(-1,-1)=Ca-1a-1

real(8),INTENT(IN) :: coeff, da0(0:3), da1(0:3)
integer,INTENT(IN) :: i,j,k
real(8) :: Caa(-2:0,-2:0), Ci(3), Cj(3), Ck(3)
real(8) :: fij(3), fik(3), fjk(3), rij(3), rik(3), rjk(3)
real(8) :: CCisqr, coCC
real(8) :: dr(3), ff(3)

Caa( 0,0) = da0(0)**2; Caa( 0,-1) = sum(da0(1:3)*da1(1:3))
Caa(-1,0) = Caa(0,-1); Caa(-1,-1) = da1(0)**2

CCisqr = 1.d0/( da0(0)*da1(0) )
coCC = coeff*CCisqr

!--- Some of calculations are unnecessary due to the action-reaction relation.
Ci(1) = -( Caa(0,-1)/Caa(0,0) )
Ci(2) =  1.d0
!Cj(1) =  Caa( 0,-1)/Caa( 0, 0) + 1.d0 
!Cj(2) = -( Caa(0,-1)/Caa(-1,-1) + 1.d0 ) 
Ck(1) = -1.d0
Ck(2) =  Caa(0,-1)/Caa(-1,-1)

rij(1:3) = da0(1:3)
rjk(1:3) = da1(1:3)

fij(1:3) = coCC*(Ci(1)*rij(1:3) + Ci(2)*rjk(1:3))
fjk(1:3) =-coCC*(Ck(1)*rij(1:3) + Ck(2)*rjk(1:3))

f(1:3,i) = f(1:3,i) + fij(1:3)
f(1:3,j) = f(1:3,j) - fij(1:3) + fjk(1:3)
f(1:3,k) = f(1:3,k) - fjk(1:3)

#ifdef STRESS
ia=i; ja=j; dr(1:3)=rij(1:3); ff(1:3)=-fij(1:3)
include 'stress'
ia=j; ja=k; dr(1:3)=rjk(1:3); ff(1:3)=-fjk(1:3)
include 'stress'
#endif

!--- Check N3rd ---
!print'(a,6f20.13)','N3rd: ', &
!Ci(1)-Ci(2), -Cj(1), Ci(2), -Ck(1), Cj(2), Ck(1)-Ck(2)

end subroutine

!-----------------------------------------------------------------------
subroutine cross_product(dr1, dr2, crs)
use atoms
! Calculate a cross product <dr1(1:3)> x <dr2(1:3)> = <crs(1:3)>
! <dr1> and <dr2> must have thier norm in 0th element.
!-----------------------------------------------------------------------
   implicit none
   real(8),INTENT(IN) :: dr1(0:3), dr2(0:3)
   real(8),INTENT(OUT) :: crs(0:3)
   real(8)  :: ndr1(1:3), ndr2(1:3)

   ndr1(1:3) = dr1(1:3)/dr1(0)
   ndr2(1:3) = dr2(1:3)/dr2(0)

   crs(1) = ndr1(2)*ndr2(3) - ndr1(3)*ndr2(2)
   crs(2) = ndr1(3)*ndr2(1) - ndr1(1)*ndr2(3)
   crs(3) = ndr1(1)*ndr2(2) - ndr1(2)*ndr2(1)
   crs(0) = sqrt( sum(crs(1:3)*crs(1:3)) )
   if(crs(0)<NSMALL) crs(0) = NSMALL
   
end subroutine 

END subroutine FORCE
