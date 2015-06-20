!-------------------------------------------------------------------------------------------
SUBROUTINE GETPARAMS()
use parameters; use atoms
!-------------------------------------------------------------------------------------------
!  This subroutine is designed solely to obtain the parameters used in the Ecalc.f90 
!    program from the rxmda.in input file.  It is similar to the input routine used
!    by Adri in "reac.f::ffinput"
!-------------------------------------------------------------------------------------------
implicit none

integer :: i,j,k,l,n,m,p,inxn   !counters for initialization
integer :: i0,i1,i2,i3,i4,ih  !Counters: # corresp to #-atom depend 

!--- Parameters that count number of entries in each field:
!integer :: nso     !* of atom types (in mod) 
!integer :: nboty  !# of bond terms given (in mod)
integer :: nodmty  !# of off-diag terms given   
integer :: npar, nvaty, ntoty, nhbty

!--- Readin Fields converted to Parameters:
integer :: nodm1, nodm2  !ID bonds to alter in off-diag terms
real(8) :: deodmh,rodmh,godmh ! off-diag term of Evdw parameters
real(8) :: rsig,rpi,rpi2 !temp storage for r0s, etc terms in off-diag part
integer :: nso2          !nso (# of types of atoms) squared 
integer :: typea,typeb   !Temp storage for filling inxn2(:,:) table

!--- NULL Transfer Fields not needed in program (used to calc other values): 
real(8) :: dnull
real(8),allocatable :: vpar(:), bo131(:), bo132(:), bo133(:)
real(8),allocatable :: rvdw1(:), eps(:), alf(:), vop(:)

dnull = 0.d0
!--- Start Getting Parameters
open(4,file=trim(adjustl(FFPath)),status="old")
read(4,'(a70)') pfile

read(4,*) npar  !num of parameters (independ of atom choice)

allocate(vpar(npar),stat=ast)

do i0=1, npar 
   read(4,1300) vpar(i0)  !temp variable...some apparently depend on atype
enddo
  
!--- Constant parameters reset to actual use:
pvdW1 = vpar(29)
pvdW1h = 0.5d0*pvdW1 
pvdW1inv = 1.d0/pvdW1

!--- a small modification in sigma-bond prime <kn>
vpar30 = vpar(30)  

!--- Parameters with 1-atom Depend,  Nr of types of atoms
read(4,'(i3)') nso    

!--- Allocation of variables:
allocate(rat(nso),rapt(nso),vnq(nso),stat=ast)
allocate(r0s(nso,nso),r0p(nso,nso),r0pp(nso,nso), stat=ast)
allocate(Val(nso),Valboc(nso),stat=ast)
allocate(mass(nso),stat=ast)
allocate(bo131(nso),bo132(nso),bo133(nso),stat=ast)
allocate(inxn2(nso,nso),inxn3(nso,nso,nso),inxn3hb(nso,nso,nso),inxn4(nso,nso,nso,nso),stat=ast)
allocate(atmname(nso), stat=ast)
allocate(Vale(nso), plp1(nso), nlpopt(nso), plp2(nso), stat=ast)
allocate(povun2(nso),povun3(nso),povun4(nso),stat=ast)
allocate(povun5(nso),povun6(nso),povun7(nso),povun8(nso),stat=ast)
!--- Valency Terms (j-dependancy only):
allocate(pval3(nso),pval5(nso), Valangle(nso),Valval(nso),stat=ast)
!--- Van der Waals Terms:
allocate(rvdw1(nso), rvdW(nso, nso), eps(nso), Dij(nso,nso), stat=ast)
allocate(alf(nso), alpij(nso,nso), stat=ast)
allocate(vop(nso), gamW(nso,nso), stat=ast)

!--- Coulomb & Charge equilibration:
allocate(chi(nso), eta(nso), gam(nso), gamij(nso,nso), stat=ast)
 
!--- Parameters that still don't depend on atom type yet
plp1(1:nso) = vpar(16)
povun3(1:nso) = vpar(33)
povun4(1:nso) = vpar(32)
povun6(1:nso) = vpar(7)
povun7(1:nso) = vpar(9)
povun8(1:nso) = vpar(10)

!--- skip 3 lines
read(4,*)
read(4,*)
read(4,*)

do i1=1, nso  !collect info on each type of atom
   read(4,1200) atmname(i1),rat(i1),Val(i1),mass(i1),rvdw1(i1),eps(i1),gam(i1),rapt(i1),Vale(i1)
   read(4,1250) alf(i1),vop(i1),Valboc(i1),povun5(i1),dnull,chi(i1),eta(i1),dnull
   read(4,1250) vnq(i1),plp2(i1),dnull,bo131(i1),bo132(i1),bo133(i1),dnull,dnull   
   read(4,1250) povun2(i1), pval3(i1),dnull,Valval(i1),pval5(i1)
enddo

nlpopt(1:nso) = 0.5d0*(Vale(1:nso) - Val(1:nso))
!--- duplicate values
Valangle(1:nso) = Valboc(1:nso)

!--- Calc default r0s, r0p, r0pp:
do i1=1,nso
   do i2=1,nso
!--- Terms for the Bond Order Calculation:
      r0s(i1,i2) = 0.5d0*(rat(i1)+rat(i2))
      r0p(i1,i2) = 0.5d0*(rapt(i1)+rapt(i2))
      r0pp(i1,i2) = 0.5d0*(vnq(i1)+vnq(i2))   
    
!--- Terms used in van der Waals calc: 
      rvdW(i1,i2) = sqrt( 4.d0*rvdw1(i1)*rvdw1(i2) )
      Dij(i1,i2) = sqrt(eps(i1)*eps(i2))
      alpij(i1,i2) = sqrt( alf(i1)*alf(i2) )
      gamW(i1,i2) = sqrt( vop(i1)*vop(i2) )  
      gamij(i1,i2) = ( gam(i1)*gam(i2) )**(-1.5d0) !<- gamcco in reac.f
   enddo
enddo  

!--- Parameters with 2-atom Depend:
read(4,1100) nboty  !# of bonds' params given 

!--- Allocation of variables:
allocate(pbo1(nboty),pbo2(nboty),pbo3(nboty),pbo4(nboty),pbo5(nboty),pbo6(nboty),bom(nboty),stat=ast)
allocate(pboc1(nboty),pboc2(nboty),pboc3(nboty),pboc4(nboty),pboc5(nboty),stat=ast)
allocate(desig(nboty), depi(nboty),depipi(nboty),pbe1(nboty),pbe2(nboty),stat=ast) 
allocate(povun1(nboty),ovc(nboty), v13cor(nboty),stat=ast)

!--- skip one line
read(4,*)

inxn2(1:nso,1:nso) = 0   !allows later flag to tell when a combination is not allowed
ih=0
do i2=1,nboty
    ih=ih+1
    read(4,1400) typea,typeb,Desig(ih),Depi(ih),Depipi(ih),pbe1(ih),pbo5(ih),v13cor(ih),pbo6(ih),povun1(ih)
    read(4,1450) pbe2(ih),pbo3(ih),pbo4(ih),bom(ih),pbo1(ih),pbo2(ih),ovc(ih),dnull  
    inxn2(typea,typeb) = ih
    inxn2(typeb,typea) = ih
enddo


!--- TEMP (required by input file backwards setup)
pboc1(1:nboty) = vpar(1)
pboc2(1:nboty) = vpar(2)

!--- for debugging <kn>
vpar1 = vpar(1)
vpar2 = vpar(2)
 
do i1=1,nso
  do i2=1,nso 
    inxn = inxn2(i1,i2)
    if(inxn/=0) then 
       pboc3(inxn) = sqrt(bo132(i1)*bo132(i2)) ! be careful about variable name
       pboc4(inxn) = sqrt(bo131(i1)*bo131(i2)) ! bo132 -> pboc2, bo131->pbo4  
       pboc5(inxn) = sqrt(bo133(i1)*bo133(i2))
      endif
  enddo
enddo


!--- Changes to off-diagonal terms:
read(4,1100) nodmty  !# of off-diag terms 
do i2=1, nodmty
   read(4,1400) nodm1,nodm2,deodmh,rodmh,godmh,rsig,rpi,rpi2
   if(rsig.GT.0.d0) r0s(nodm1,nodm2)=rsig
   if(rsig.GT.0.d0) r0s(nodm2,nodm1)=rsig 
   if(rpi.GT.0.d0)  r0p(nodm1,nodm2)=rpi
   if(rpi.GT.0.d0)  r0p(nodm2,nodm1)=rpi
   if(rpi2.GT.0.d0) r0pp(nodm1,nodm2)=rpi2
   if(rpi2.GT.0.d0) r0pp(nodm2,nodm1)=rpi2
   if (rodmh.GT.0.d0) rvdW(nodm1,nodm2)=2.0*rodmh 
   if (rodmh.GT.0.d0) rvdW(nodm2,nodm1)=2.0*rodmh
   if (deodmh.GT.0.d0) Dij(nodm1,nodm2)=deodmh
   if (deodmh.GT.0.d0) Dij(nodm2,nodm1)=deodmh
   if (godmh.GT.0.d0) alpij(nodm1,nodm2)=godmh
   if (godmh.GT.0.d0) alpij(nodm2,nodm1)=godmh
enddo

!--- Input Valency Terms from Input File
inxn3(1:nso,1:nso,1:nso) = 0
read(4,1100) nvaty

allocate(pval1(nvaty),pval2(nvaty),pval4(nvaty),pval6(nvaty),stat=ast)
allocate(pval7(nvaty),pval8(nvaty),pval9(nvaty),pval10(nvaty),stat=ast)
allocate(theta00(nvaty),stat=ast)
allocate(ppen1(nvaty),ppen2(nvaty),ppen3(nvaty),ppen4(nvaty),stat=ast)
allocate(pcoa1(nvaty),pcoa2(nvaty),pcoa3(nvaty),pcoa4(nvaty),stat=ast)

do i=1, nvaty
   read(4,1500) i1,i2,i3,theta00(i),pval1(i),pval2(i),pcoa1(i),pval7(i),ppen1(i),pval4(i)
   inxn3(i1,i2,i3) = i
   inxn3(i3,i2,i1) = i
enddo

inxn3(1,2,2)=0 !react.f, line 4933
inxn3(2,2,1)=0 !react.f, line 4933

!--- Valency Terms which do not depend on inxn type:
pval6(1:nvaty) = vpar(15)
pval8(1:nvaty) = vpar(34)
pval9(1:nvaty) = vpar(17)
pval10(1:nvaty) = vpar(18)
!--- Penalty Terms which do not depend on inxn type:
ppen2(1:nvaty) = vpar(20)
ppen3(1:nvaty) = vpar(21)
ppen4(1:nvaty) = vpar(22)
!--- 3body Conjugation Terms which do not depend on type:
pcoa2(1:nvaty) = vpar(3)
pcoa3(1:nvaty) = vpar(39)
pcoa4(1:nvaty) = vpar(31)
!--- theta00 given in degrees, but used in radians. Convert by:
theta00(1:nvaty) = (pi/180.d0)*theta00(1:nvaty) 


read(4,1100) ntoty
allocate(ptor1(ntoty),ptor2(ntoty),ptor3(ntoty),ptor4(ntoty),V1(ntoty), V2(ntoty),V3(ntoty),stat=ast)
allocate(pcot1(ntoty),pcot2(ntoty),stat=ast)

inxn4(1:nso,1:nso,1:nso,1:nso) = 0  
do i=1,ntoty
   read(4,1600)i1,i2,i3,i4,V1(i),V2(i),V3(i),ptor1(i),pcot1(i),dnull,dnull 
!--- Set up inxn4 lookup reference array
      if(i1==0) then   !condensed input, means that all i1,i4 for this arrangement of i2,i3 are the same
         do i1=1,nso
         do i4=1,nso
         if(inxn4(i1,i2,i3,i4)==0.and.inxn4(i1,i3,i2,i4)==0) then
            inxn4(i1,i2,i3,i4) = i
            inxn4(i4,i2,i3,i1) = i
            inxn4(i1,i3,i2,i4) = i
            inxn4(i4,i3,i2,i1) = i
          endif
          enddo
          enddo
      else
          inxn4(i1,i2,i3,i4) = i
          inxn4(i4,i2,i3,i1) = i
          inxn4(i1,i3,i2,i4) = i
          inxn4(i4,i3,i2,i1) = i 
      endif
enddo

!and a few which don't depend on type 
ptor2(1:ntoty) = vpar(24)
ptor3(1:ntoty) = vpar(25)
ptor4(1:ntoty) = vpar(26)  
pcot2(1:ntoty) = vpar(28)

!--- Input Hydrogen Bond Terms
inxn3hb(1:nso,1:nso,1:nso) = 0
read(4,1100) nhbty
allocate(phb1(nhbty),phb2(nhbty),phb3(nhbty),r0hb(nhbty),stat=ast); if(ast/=0) stop "Error: H-bond Terms@init"

do i=1,nhbty
   read(4,1500) i1,i2,i3,r0hb(i),phb1(i),phb2(i),phb3(i)
   inxn3hb(i1,i2,i3) = i    !Note: inxn3hb(i,j,k) /= inxn3hb(k,j,i)
enddo

!--- close parameter file "ffield"
close(4)

!--- Formats:
1100 format (i3,2x,a2,3x,3d22.15)
1200 format (1x,a2,10f9.4)
1250 format (3x,10f9.4)
1300 format (f10.4)
1400 format (2i3,8f9.4)
1450 format (6x,8f9.4)
1500 format (3i3,7f9.4)
1600 format (4i3,7f9.4)


!--- coefficient of coulomb energy
Cclmb = Cclmb0  !Eclmb

!--- Parameters for charge variable routine. 
!--- In original parameter, chiEEM and etaEEM are given in [ev], not [kcal/mol]
!--- Definition of the stiffness parameter <eta> is different from 
!--- the original code and our code. It's need to be multiplied by 2.
eta(:) = eta(:)*2.d0

rctap = rctap0
rctap2 = rctap**2

CTap(0:7)=(/1.d0, 0.d0, 0.d0, 0.d0,   -35.d0/(rctap)**4, &
          84.d0/(rctap)**5, -70.d0/(rctap)**6, &
          20.d0/(rctap)**7 /)

END SUBROUTINE
