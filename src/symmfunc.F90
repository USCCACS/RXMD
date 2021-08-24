module symmfunc

  implicit none

  type sf_F1_return 
    double precision :: Fijk
    double precision, dimension(3) :: dFijk_dRj
    double precision, dimension(3) :: dFijk_dRk
  end type sf_F1_return

  double precision, parameter, private :: PI    = 3.14159265358979d0
  double precision, parameter, private :: PI2   = 2.0d0*PI
  double precision, parameter, private :: EPS   = 1.0d-12
  integer,          parameter, private :: NG    = 5

  !--------------------------------------------------------------------!
  ! sf_nG_type(i)    : number of symmetry functions for species i      !
  !                                                                    !
  ! sf_nTypes        : number of different atom/point types            !
  ! sf_nPairs        : number type pairs (= nTypes!)                   !
  ! sf_nGmax         : max. num. symmetry functs. per type pair/triple !
  ! sf_idx2(i,j)     : pair index for species i-j                      !
  ! sf_idx3(i,j,k)   : triple index for species i-j-k                  !
  ! sf_nG(i,pij)     : num. symmetry function type i for type pair pij !
  ! sf_pG?(i,iG,pij) : parameter i for the iG-th symm.-func. of type ? !
  !                    for type pair pij                               !
  !                    Parameters:  G1:  Rc                            !
  !                                 G2:  Rc, Rs, eta                   !
  !                                 G3:  Rc, kappa                     !
  !                                 G4:  Rc, lambda, zeta, eta         !
  !                                 G5:  Rc, lambda, zeta, eta         !
  !--------------------------------------------------------------------!

  integer,          dimension(:),     allocatable, public  :: sf_nG_type

  integer,                                         private :: sf_nTypes
  integer,                                         private :: sf_nPairs
  integer,                                         private :: sf_nTriples
  integer,                                         private :: sf_nGMax
  integer,          dimension(:,:),   allocatable, private :: sf_idx2
  integer,          dimension(:,:,:), allocatable, private :: sf_idx3
  integer,          dimension(:),     allocatable, private :: sf_iG02
  integer,          dimension(:),     allocatable, private :: sf_iG03
  integer,          dimension(:,:),   allocatable, private :: sf_nG2
  integer,          dimension(:,:),   allocatable, private :: sf_nG3
  double precision, dimension(:,:),   allocatable, private :: sf_pG1
  double precision, dimension(:,:,:), allocatable, private :: sf_pG2
  double precision, dimension(:,:,:), allocatable, private :: sf_pG3
  double precision, dimension(:,:,:), allocatable, private :: sf_pG4
  double precision, dimension(:,:,:), allocatable, private :: sf_pG5

  !$omp declare target (sf_pG4)
  !--------------------------------------------------------------------!

  integer, public :: stdout = 6
  integer, public :: stderr = 0

  !--------------------------------------------------------------------!

  logical, private :: isInit  = .false.
  integer, private :: memSize = 0

  !--------------------------------------------------------------------!

  integer, parameter :: norder = 1
  integer :: nbin = 100000 ! size of table
!  integer, parameter :: nbin = 10000 ! for the second-order interpolation
  real(8), allocatable, private :: cutv(:), cutva(:), cutd(:), cutda(:)
  real(8) :: dx, dxr

  real(8), allocatable, private :: exv(:), exva(:)
  real(8) :: exx, dex, dexr
  !$omp declare target (exx, dexr, exv)
  save

contains !-------------------------------------------------------------!

  subroutine sf_init(ntypes, nGmax)

    implicit none

    integer, intent(in) :: ntypes
    integer, intent(in) :: nGmax

    integer :: i, j, k, idx2, idx3

    sf_nTypes = ntypes
    sf_nGMax  = nGmax

    ! initialize type pair/triple index:
    ! note: for angular symmetry functions A-B-C = A-C-B
    allocate(sf_idx2(ntypes,ntypes), sf_idx3(ntypes,ntypes,ntypes))
    idx2 = 0
    idx3 = 0
    do i = 1, ntypes
       do j = 1, ntypes
          idx2 = idx2 + 1
          sf_idx2(i,j) = idx2
          idx3 = idx3 + 1
          sf_idx3(i,j,j) = idx3
          do k = j+1, ntypes
             idx3 = idx3 + 1
             sf_idx3(i,j,k) = idx3
             sf_idx3(i,k,j) = idx3
          end do
       end do
    end do
    sf_nPairs = idx2
    sf_nTriples = idx3

    allocate(sf_nG_type(ntypes),            &
             sf_nG2(NG, sf_nPairs),         &
             sf_nG3(NG, sf_nTriples),       &
             sf_iG02(sf_nPairs),            &
             sf_iG03(sf_nTriples),          &
             sf_pG1(nGmax, sf_nPairs),      &
             sf_pG2(3, nGmax, sf_nPairs),   &
             sf_pG3(2, nGmax, sf_nPairs),   &
             sf_pG4(4, nGmax, sf_nTriples), &
             sf_pG5(4, nGmax, sf_nTriples)  )
    sf_nG_type(:) = 0
    sf_nG2(:,:) = 0
    sf_nG3(:,:) = 0

    isInit = .true.

  end subroutine sf_init

  !--------------------------------------------------------------------!
  !                      parameter initialization                      !
  !--------------------------------------------------------------------!

  subroutine sf_add_rad(funct, type1, type2, Rc, Rs, eta, kappa)

    implicit none

    integer,                    intent(in) :: funct
    integer,                    intent(in) :: type1, type2
    double precision, optional, intent(in) :: Rc, Rs, eta, kappa

    integer :: iG, pij

    pij = sf_idx2(type1,type2)

    if (sf_nG2(funct,pij) >= sf_nGMax) then
       write(stderr,*) "Error: Memory exceeded in `sf_add()'."
       write(stderr,*) sf_nG2(funct,pij), sf_nGMax
       stop
    end if

    if (.not. present(Rc)) then
       write(stderr,*) "Error: Rc undefined in `sf_add()'"
       stop
    end if

    iG = sf_nG2(funct,pij) + 1
    select case(funct)
    case(1)
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG1(iG,pij)   = Rc
    case(2)
       if (.not. (present(Rs) .and. present(eta))) then
          write(stderr,*) "Error: needed for G2: Rs, eta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG2(1,iG,pij) = Rc
       sf_pG2(2,iG,pij) = Rs
       sf_pG2(3,iG,pij) = eta
    case(3)
       if (.not. present(kappa)) then
          write(stderr,*) "Error: needed for G3: kappa"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG3(1,iG,pij) = Rc
       sf_pG3(2,iG,pij) = kappa
    case default
       write(stderr,*) "Error: Not a radial symmetry function type: ", funct
       stop
    end select

    call sf_build_index()

  end subroutine sf_add_rad

  !--------------------------------------------------------------------!

  subroutine sf_add_ang(funct, type1, type2, type3, Rc, eta, lambda, zeta)

    implicit none

    integer,                    intent(in) :: funct
    integer,                    intent(in) :: type1, type2, type3
    double precision, optional, intent(in) :: Rc, eta, lambda, zeta

    integer :: iG, tijk

    tijk = sf_idx3(type1, type2, type3)

    if (sf_nG3(funct,tijk) >= sf_nGMax) then
       write(stderr,*) "Error: Memory exceeded in `sf_add()'."
       write(stderr,*) sf_nG3(funct,tijk), sf_nGMax
       stop
    end if

    if (.not. present(Rc)) then
       write(stderr,*) "Error: Rc undefined in `sf_add()'"
       stop
    end if

    iG = sf_nG3(funct,tijk) + 1
    select case(funct)
    case(4)
       if (.not. (present(lambda) .and. present(zeta) &
           .and. present(eta))) then
          write(stderr,*) "Error: needed for G4: eta, lambda, zeta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG3(funct,tijk) = iG
       sf_pG4(1,iG,tijk) = Rc
       sf_pG4(2,iG,tijk) = lambda
       sf_pG4(3,iG,tijk) = zeta
       sf_pG4(4,iG,tijk) = eta
    case(5)
       if (.not. (present(lambda) .and. present(zeta) &
           .and. present(eta))) then
          write(stderr,*) "Error: needed for G5: eta, lambda, zeta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG3(funct,tijk) = iG
       sf_pG5(1,iG,tijk) = Rc
       sf_pG5(2,iG,tijk) = lambda
       sf_pG5(3,iG,tijk) = zeta
       sf_pG5(4,iG,tijk) = eta
    case default
       write(stderr,*) "Error: Not an angular symmetry function type: ", funct
       stop
    end select

    call sf_build_index()

  end subroutine sf_add_ang

  !--------------------------------------------------------------------!
  ! The functions are organized in memory not in the same order as     !
  ! they are added, but rather depending on their type and species.    !
  ! The procedure `sf_build_index()' generates the indices for radial  !
  ! and angular functions.                                             !
  !--------------------------------------------------------------------!

  subroutine sf_build_index()

    implicit none

    integer :: pij, tijk, i
    integer :: itype1, itype2, itype3

    do itype1 = 1, sf_nTypes
       i = 0
       do itype2 = 1, sf_nTypes
          pij = sf_idx2(itype1, itype2)
          sf_iG02(pij) = i
          i = i + sum(sf_nG2(:,pij))
          do itype3 = itype2, sf_nTypes
             tijk = sf_idx3(itype1, itype2, itype3)
             sf_iG03(tijk) = i
             i = i + sum(sf_nG3(:,tijk))
          end do
       end do
    end do

  end subroutine sf_build_index

  !--------------------------------------------------------------------!
  !               evaluate finger print for single atom                !
  !--------------------------------------------------------------------!

  subroutine sf_fingerprint(type1, Xi, nx, X, typej, n, G, dGi, dGj, dGh)

    implicit none

    !------------------------------------------------------------------!
    ! type1     : type index of the central particle                   !
    ! Xi(1:3)   : cartesian coordinates of the central particle        !
    ! nx        : number of other particles                            !
    ! X(1:3,j)  : cartesian coordinates of the j-th particle           !
    ! typej(j)  : atom type index of j-th atom                         !
    ! n         : number of symmetry functions                         !
    ! G(j)      : (out) value of the j-th symmetry function            !
    ! dGi(1:3,j): (out, optional) derivative of the j-th symmetry      !
    !             function wrt. the central particle's coordinates     !
    ! dGj(:,j,k): (out, optional) derivative of the j-th symmetry      !
    !             function wrt. the coordinates of the k-th particle   !
    ! Note: dGi and dGj have to be present/absent simultaneously       !
    !------------------------------------------------------------------!

    integer,                                       intent(in)  :: type1
    double precision, dimension(3),                intent(in)  :: Xi
    integer,                                       intent(in)  :: nx
    double precision, dimension(3,nx),             intent(in)  :: X
    integer,          dimension(nx),               intent(in)  :: typej
    integer,                                       intent(in)  :: n
    double precision, dimension(n),                intent(out) :: G
    double precision, dimension(3,n),    optional, intent(out) :: dGi
    double precision, dimension(3,n,nx), optional, intent(out) :: dGj
    double precision, dimension(6,n),    optional, intent(out) :: dGh

    integer          :: j, k, pij, tijk
    integer          :: iG_ang0, iG
    double precision :: Rij, Rik, Rjk, cost
    logical          :: deriv, lstrs

    double precision, dimension(3) :: R1, R2, R3

    lstrs = .false.
    if (present(dGi) .and. present(dGj)) then
       deriv = .true.
       if (present(dGh)) then
           lstrs = .true.
       end if
    else
       deriv = .false.
    end if

    if (.not. isInit) return

    if (n < sf_nG_type(type1)) then
       write(stderr,*) "Error: number of symmetry functions exceeded."
       stop
    end if

    G(1:n)  = 0.0d0
    iG_ang0 = 0
    if (deriv) then
       dGi(1:3,1:n) = 0.0d0
       dGj(1:3,1:n,1:nx) = 0.0d0
       if (lstrs) then
          dGh(1:6,1:n) = 0.0d0
       end if
    end if

    !$omp parallel do
    jloop : do j = 1, nx
       pij = sf_idx2(type1, typej(j))

       R1  = X(:,j) - Xi(:)
       Rij = sqrt(sum(R1*R1))
       R1  = R1/Rij

       iG = sf_iG02(pij)

       ! radial symmetry functions:
       if (deriv) then
          if (lstrs) then
             call sf_G1_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j), dGh)
             call sf_G2_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j), dGh)
             call sf_G3_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j), dGh)
          else
             call sf_G1_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
             call sf_G2_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
             call sf_G3_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
          end if
       else
          call sf_G1_update(pij, R1(1:3), Rij, n, G, iG)
          call sf_G2_update(pij, R1(1:3), Rij, n, G, iG)
          call sf_G3_update(pij, R1(1:3), Rij, n, G, iG)
       end if

       ! if (iG_ang0 == 0) iG_ang0 = sum(sf_nG2(:,:))
       kloop : do k = j+1, nx
          tijk = sf_idx3(type1, typej(j), typej(k))

          R2  = X(:,k) - Xi(:)
          Rik = sqrt(sum(R2*R2))
          if (Rik <= 1.0d-8) then
             write(stderr, *) "Warning: an atom in the neighbor list is " &
                  // "identical to the central atom (sf_fingerprint)"
          else
             R2  = R2/Rik
          end if

          R3  = X(:,k) - X(:,j)
          Rjk = sqrt(sum(R3*R3))
          if (Rjk <= 1.0d-8) then
             write(stderr, *) "Warning: redundant atoms in neighbor " &
                              // "list (sf_fingerprint)"
          else
             R3  = R3/Rjk
          end if

          ! cos(theta_ijk)
          cost = sum(R1*R2)
          cost = max(cost,-1.0d0)
          cost = min(cost, 1.0d0)

          iG = sf_iG03(tijk)

          ! angular symmetry functions:
          if (deriv) then
             if (lstrs) then
                call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k), dGh)
                call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k), dGh)
             else
                call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k))
                call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k))
             end if
          else
             call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, &
                               cost, n, G, iG)
             call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                               n, G, iG)
          end if

       end do kloop
    end do jloop

  end subroutine sf_fingerprint



  !==================== radial symmetry functions =====================!



  !--------------------------------------------------------------------!
  !            radial symm. function 1 (eq. 5 in Ref. [1])             !
  !                                                                    !
  ! here: vecRij = vec{R}_ij/||vec{R})_ij||                            !
  !--------------------------------------------------------------------!

  subroutine sf_G1_ij(Rij, Rc, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc, Rcr

    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

    Gij  = fc
    dGij = dfc

  end subroutine sf_G1_ij

  !--------------------------------------------------------------------!

  subroutine sf_G1_update(pij, vecRij, Rij, n, G, iG, dGi, dGj, dGh)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer          :: iG1
    double precision :: G1, dG1, dG1r
    double precision :: Rc

    do iG1 = 1, sf_nG2(1,pij)
       iG    = iG + 1
       Rc    = sf_pG1(iG1,pij)
       call sf_G1_ij(Rij, Rc, G1, dG1)
       G(iG) = G(iG) + G1
       if (present(dGi)) dGi(1:3,iG) = dGi(1:3,iG) - dG1*vecRij(1:3)
       if (present(dGj)) dGj(1:3,iG) = dGj(1:3,iG) + dG1*vecRij(1:3)
       if (present(dGh)) then
           dG1r = dG1*Rij
           dGh(1,iG) = dGh(1,iG) + dG1r*vecRij(1)*vecRij(1)
           dGh(2,iG) = dGh(2,iG) + dG1r*vecRij(2)*vecRij(2)
           dGh(3,iG) = dGh(3,iG) + dG1r*vecRij(3)*vecRij(3)
           dGh(4,iG) = dGh(4,iG) + dG1r*vecRij(2)*vecRij(3)
           dGh(5,iG) = dGh(5,iG) + dG1r*vecRij(3)*vecRij(1)
           dGh(6,iG) = dGh(6,iG) + dG1r*vecRij(1)*vecRij(2)
       end if
    end do

  end subroutine sf_G1_update

  !--------------------------------------------------------------------!
  !            radial symm. function 2 (eq. 6 in Ref. [1])             !
  !--------------------------------------------------------------------!

  subroutine sf_G2_ij(Rij, Rc, Rs, eta, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(in)  :: Rs
    double precision, intent(in)  :: eta
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc, Rcr
    double precision :: fexp, arg, x, d, etag
    integer :: m

    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

    arg  = Rij - Rs
!    fexp = exp(-eta*arg*arg)
    etag = eta*arg
    x = etag*arg
    if( x > exx ) then
        fexp = 0.d0
    else
    if( norder == 1 ) then
        !---first-order interpolation
        d = x*dexr
        m = d
        d = d - m
        fexp = (1d0-d)*exv(m)+d*exv(m+1)
    else
        !---second-order interpolation
        m = x/(2.d0*dex)
        m = 2*m
        d = 0.5d0*( x*dexr - dble(m) )
        fexp = d*( (d-1.d0)*exva(m) + exv(m+2) - exv(m) ) + exv(m)
    end if
    end if

    Gij  = fexp*fc
    dGij = fexp*( dfc - 2.0d0*etag*fc  )

  end subroutine sf_G2_ij

  !--------------------------------------------------------------------!

  subroutine sf_G2_update(pij, vecRij, Rij, n, G, iG, dGi, dGj, dGh)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer          :: iG2
    double precision :: Rc, Rs, eta
    double precision :: G2, dG2, vecRij2(3)

    do iG2 = 1, sf_nG2(2,pij)
       iG    = iG + 1
       Rc    = sf_pG2(1,iG2,pij)
       Rs    = sf_pG2(2,iG2,pij)
       eta   = sf_pG2(3,iG2,pij)
       call sf_G2_ij(Rij, Rc, Rs, eta, G2, dG2)
       G(iG) = G(iG) + G2
       if (present(dGi)) then
           vecRij2(1:3) = dG2*vecRij(1:3)
           dGi(1:3,iG) = dGi(1:3,iG) - vecRij2(1:3)
!       end if
!       if (present(dGj)) then
           dGj(1:3,iG) = dGj(1:3,iG) + vecRij2(1:3)
!       end if
         if (present(dGh)) then
             vecRij2(1:3) = vecRij2(1:3)*Rij
             dGh(1,iG) = dGh(1,iG) + vecRij2(1)*vecRij(1)
             dGh(2,iG) = dGh(2,iG) + vecRij2(2)*vecRij(2)
             dGh(3,iG) = dGh(3,iG) + vecRij2(3)*vecRij(3)
             dGh(4,iG) = dGh(4,iG) + vecRij2(2)*vecRij(3)
             dGh(5,iG) = dGh(5,iG) + vecRij2(3)*vecRij(1)
             dGh(6,iG) = dGh(6,iG) + vecRij2(1)*vecRij(2)
         end if
       end if
    end do

  end subroutine sf_G2_update

  !--------------------------------------------------------------------!
  !            radial symm. function 3 (eq. 6 in Ref. [1])             !
  !--------------------------------------------------------------------!

  subroutine sf_G3_ij(Rij, Rc, kappa, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(in)  :: kappa
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc, Rcr
    double precision :: fcos, fsin, arg

    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

    arg  = kappa*Rij
    fcos = cos(arg)
    fsin = sin(arg)

    Gij  = fcos*fc
    dGij = fcos*dfc - kappa*fsin*fc

  end subroutine sf_G3_ij

  !--------------------------------------------------------------------!

  subroutine sf_G3_update(pij, vecRij, Rij, n, G, iG, dGi, dGj, dGh)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer          :: iG3
    double precision :: Rc, kappa
    double precision :: G3, dG3, dG3r

    do iG3 = 1, sf_nG2(3,pij)
       iG    = iG + 1
       Rc    = sf_pG3(1,iG3,pij)
       kappa = sf_pG3(2,iG3,pij)
       call sf_G3_ij(Rij, Rc, kappa, G3, dG3)
       G(iG) = G(iG) + G3
       if (present(dGi)) dGi(1:3,iG) = dGi(1:3,iG) - dG3*vecRij(1:3)
       if (present(dGj)) dGj(1:3,iG) = dGj(1:3,iG) + dG3*vecRij(1:3)
       if (present(dGh)) then
           dG3r = dG3*Rij
           dGh(1,iG) = dGh(1,iG) + dG3r*vecRij(1)*vecRij(1)
           dGh(2,iG) = dGh(2,iG) + dG3r*vecRij(2)*vecRij(2)
           dGh(3,iG) = dGh(3,iG) + dG3r*vecRij(3)*vecRij(3)
           dGh(4,iG) = dGh(4,iG) + dG3r*vecRij(2)*vecRij(3)
           dGh(5,iG) = dGh(5,iG) + dG3r*vecRij(3)*vecRij(1)
           dGh(6,iG) = dGh(6,iG) + dG3r*vecRij(1)*vecRij(2)
       end if
    end do

  end subroutine sf_G3_update



  !==================== angular symmetry functions ====================!

  !--------------------------------------------------------------------!
  !          first angular symm. function (eq. 8 in Ref. [1])          !
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !                        angle dependend part                        !
  ! This is function F_1(R_ij,R_ik) in the documentation.              !
  !--------------------------------------------------------------------!
  
  !$omp declare target (sf_F1_ijk_func)
  function sf_F1_ijk_func(vecRij, vecRik, Rijr, Rikr, cost, lambda, &
                       zeta)

    implicit none

    type(sf_F1_return) :: sf_F1_ijk_func
    double precision, dimension(3), intent(in)  :: vecRij, vecRik
    double precision, value,              intent(in)  :: Rijr, Rikr
    double precision, value,              intent(in)  :: cost
    double precision, value,              intent(in)  :: lambda
    double precision, value,              intent(in)  :: zeta
    !double precision,                :: Fijk
    !double precision, dimension(3),  :: dFijk_dRj
    !double precision, dimension(3),  :: dFijk_dRk

    double precision :: arg, prefactor, argz, prefactor1
    real(8), parameter :: epsi = 1.d-12

    if (abs(cost) > 1.0d0) write(*,*) "cos(theta) = ", cost

!    if( abs(lambda-1.0d0) < epsi ) then
!        arg = 0.5d0*(1.0d0 + cost)
!    else if( abs(lambda+1.0d0) < epsi ) then
!        arg = 0.5d0*(1.0d0 - cost)
!    else
        arg = 0.5d0*(1.0d0 + lambda*cost)
!    end if

!    prefactor = 0.5d0*zeta*lambda*arg**(zeta-1.0d0)
    if( abs(zeta-1.0d0) < epsi ) then
        argz = 1.d0
    else if( abs(zeta-2.0d0) < epsi ) then
        argz = arg
    else if( abs(zeta-3.0d0) < epsi ) then
        argz = arg*arg
    else if( abs(zeta-4.0d0) < epsi ) then
        argz = arg*arg*arg
    else
        argz = arg**(zeta-1.0d0)
    end if
!    if( abs(lambda-1.0d0) < epsi ) then
!        prefactor = 0.5d0*zeta*argz
!    else if( abs(lambda+1.0d0) < epsi ) then
!        prefactor = -0.5d0*zeta*argz
!    else
        prefactor = 0.5d0*zeta*lambda*argz
!    end if

!    Fijk  = arg**zeta
    sf_F1_ijk_func%Fijk  = argz*arg

    prefactor1 = prefactor*Rijr
!    dFijk_dRj(1:3) = prefactor1*( -cost*vecRij(1:3) + vecRik(1:3) )
    sf_F1_ijk_func%dFijk_dRj(1:3) = prefactor1*vecRij(1:3)

    prefactor1 = prefactor*Rikr
!    dFijk_dRk(1:3) = prefactor1*( -cost*vecRik(1:3) + vecRij(1:3) )
    sf_F1_ijk_func%dFijk_dRk(1:3) = prefactor1*vecRik(1:3)
  end function sf_F1_ijk_func


  subroutine sf_F1_ijk(vecRij, vecRik, Rijr, Rikr, cost, lambda, &
                       zeta, Fijk, dFijk_dRj, dFijk_dRk)

    implicit none

    double precision, dimension(3), intent(in)  :: vecRij, vecRik
    double precision, value,              intent(in)  :: Rijr, Rikr
    double precision, value,              intent(in)  :: cost
    double precision, value,              intent(in)  :: lambda
    double precision, value,              intent(in)  :: zeta
    double precision,               intent(out) :: Fijk
    double precision, dimension(3), intent(out) :: dFijk_dRj
    double precision, dimension(3), intent(out) :: dFijk_dRk

    double precision :: arg, prefactor, argz, prefactor1
    real(8), parameter :: epsi = 1.d-12

    if (abs(cost) > 1.0d0) write(*,*) "cos(theta) = ", cost

!    if( abs(lambda-1.0d0) < epsi ) then
!        arg = 0.5d0*(1.0d0 + cost)
!    else if( abs(lambda+1.0d0) < epsi ) then
!        arg = 0.5d0*(1.0d0 - cost)
!    else
        arg = 0.5d0*(1.0d0 + lambda*cost)
!    end if

!    prefactor = 0.5d0*zeta*lambda*arg**(zeta-1.0d0)
    if( abs(zeta-1.0d0) < epsi ) then
        argz = 1.d0
    else if( abs(zeta-2.0d0) < epsi ) then
        argz = arg
    else if( abs(zeta-3.0d0) < epsi ) then
        argz = arg*arg
    else if( abs(zeta-4.0d0) < epsi ) then
        argz = arg*arg*arg
    else
        argz = arg**(zeta-1.0d0)
    end if
!    if( abs(lambda-1.0d0) < epsi ) then
!        prefactor = 0.5d0*zeta*argz
!    else if( abs(lambda+1.0d0) < epsi ) then
!        prefactor = -0.5d0*zeta*argz
!    else
        prefactor = 0.5d0*zeta*lambda*argz
!    end if

!    Fijk  = arg**zeta
    Fijk  = argz*arg

    prefactor1 = prefactor*Rijr
!    dFijk_dRj(1:3) = prefactor1*( -cost*vecRij(1:3) + vecRik(1:3) )
    dFijk_dRj(1:3) = prefactor1*vecRij(1:3)

    prefactor1 = prefactor*Rikr
!    dFijk_dRk(1:3) = prefactor1*( -cost*vecRik(1:3) + vecRij(1:3) )
    dFijk_dRk(1:3) = prefactor1*vecRik(1:3)

  end subroutine sf_F1_ijk

  !--------------------------------------------------------------------!
  !                      distance dependend part                       !
  ! This is function F_2(R) in the documentation.                      !
  ! --> very similar to sf_G2_ij(), may be combined in future.         !
  !--------------------------------------------------------------------!
   !$omp declare target (sf_F2_ij)
   function sf_F2_ij_func(Rij, Rc, Rcr, eta)

      implicit none

      double precision, dimension(2) :: sf_F2_ij_func
      double precision,value, intent(in)  :: Rij
      double precision,value, intent(in)  :: Rc, Rcr
      double precision,value, intent(in)  :: eta
      !double precision, intent(out) :: Fij
      !double precision, intent(out) :: dFij
      
      double precision :: fc, dfc
      double precision :: fexp, x, d, etij
      integer :: m
      
   !    Rcr = 1.d0/Rc
      call sf_cut(Rij, Rc, Rcr, fc, dfc)

   !    fexp = exp(-eta*Rij*Rij)
      etij = eta*Rij
      x = etij*Rij
      if( x > exx ) then
         fexp = 0.d0
      else
      if( norder == 1 ) then
         !---first-order interpolation
         d = x*dexr
         m = d
         d = d - m
         fexp = (1d0-d)*exv(m)+d*exv(m+1)
      else
         !---second-order interpolation
         m = x/(2.d0*dex)
         m = 2*m
         d = 0.5d0*( x*dexr - dble(m) )
         fexp = d*( (d-1.d0)*exva(m) + exv(m+2) - exv(m) ) + exv(m)
      end if
      end if

      sf_F2_ij_func(1)  = fexp*fc !Fij
      sf_F2_ij_func(2) = fexp*( dfc - 2.0d0*etij*fc  ) !dFij

   end function sf_F2_ij_func

  subroutine sf_F2_ij(Rij, Rc, Rcr, eta, Fij, dFij)

    implicit none

    double precision,value, intent(in)  :: Rij
    double precision,value, intent(in)  :: Rc, Rcr
    double precision,value, intent(in)  :: eta
    double precision, intent(out) :: Fij
    double precision, intent(out) :: dFij
    
    double precision :: fc, dfc
    double precision :: fexp, x, d, etij
    integer :: m
    
!    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

!    fexp = exp(-eta*Rij*Rij)
    etij = eta*Rij
    x = etij*Rij
    if( x > exx ) then
        fexp = 0.d0
    else
      if( norder == 1 ) then
         !---first-order interpolation
         d = x*dexr
         m = d
         d = d - m
         fexp = (1d0-d)*exv(m)+d*exv(m+1)
      else
         !---second-order interpolation
         m = x/(2.d0*dex)
         m = 2*m
         d = 0.5d0*( x*dexr - dble(m) )
         fexp = d*( (d-1.d0)*exva(m) + exv(m+2) - exv(m) ) + exv(m)
      end if
    end if

    Fij  = fexp*fc
    dFij = fexp*( dfc - 2.0d0*etij*fc  )

  end subroutine sf_F2_ij

  !--------------------------------------------------------------------!

subroutine sf_G4_update(tijk, vecRij, vecRik, vecRjk, Rij, Rik, Rjk, &
                        cost, n, G, iG, dGi, dGj, dGk, dGh)

   implicit none

   integer,                                    intent(in)    :: tijk
   double precision, dimension(3),             intent(in)    :: vecRij, vecRik, vecRjk
   double precision,                           intent(in)    :: Rij, Rik, Rjk
   double precision,                           intent(in)    :: cost
   integer,                                    intent(in)    :: n
   double precision, dimension(n),             intent(inout) :: G
   integer,                                    intent(inout) :: iG
   double precision, dimension(3,n), optional, intent(inout) :: dGi
   double precision, dimension(3,n), optional, intent(inout) :: dGj
   double precision, dimension(3,n), optional, intent(inout) :: dGk
   double precision, dimension(6,n), optional, intent(inout) :: dGh

   integer                        :: iG4 
   double precision               :: Rc, Rcr, lambda, zeta, eta
   double precision               :: G4, G0
   double precision, dimension(3) :: dG4, dF1j, dF1k

   double precision :: F1, F2ij, dF2ij, F2ik, dF2ik, F2jk, dF2jk
   double precision :: Rijr, Rikr, vecRijik(3), vecRikij(3)
   double precision :: dGhh(6), dGhij(3), dGhik(3), dGhjk(3)
   double precision :: dG401, dG402, dG403, vecRij0(3), vecRik0(3), vecRjk0(3)

   type(sf_F1_return) :: sf_F1_result
   double precision, dimension(2) :: sf_F2_result

   integer :: iG_tmp
   integer :: iter 
   iter = sf_nG3(4,tijk)

    Rijr = 1.d0/Rij
    Rikr = 1.d0/Rik
    vecRijik(1:3) = -cost*vecRij(1:3) + vecRik(1:3)
    vecRikij(1:3) = -cost*vecRik(1:3) + vecRij(1:3)

    !iG = iG + iter

    !PRINT *, iter

    !omp target data map(to:sf_pG4(1:4,1:60,1:6), exx, dexr, exv(0:nbin)) 
    !$omp target teams distribute parallel do default(none) &
    !$omp&       private(iG_tmp, Rc, Rcr, lambda, zeta, eta, G0, G4, vecRij0, vecRik0, vecRjk0, &
    !$omp&               F1, dF1j, dF1k, F2ij, dF2ij, F2ik, dF2ik, F2jk, dF2jk, dG4, & 
    !$omp&               dG401, dG402, dG403, dGhij, dGhik, dGhjk, dGhh, sf_F1_result, sf_F2_result) &
    !$omp&       shared(iter,sf_pG4, tijk, G, iG, vecRijik, vecRikij, Rijr, Rikr, cost, &
    !$omp&              Rij, Rik, Rjk, dGi,dGj,dGk, vecRij, vecRik, vecRjk, dGh) 
    do iG4 = 1, iter
      iG_tmp = iG + iG4
      Rc     = sf_pG4(1,iG4,tijk)
      lambda = sf_pG4(2,iG4,tijk)
      zeta   = sf_pG4(3,iG4,tijk)
      eta    = sf_pG4(4,iG4,tijk)

      !call sf_F1_ijk(vecRijik, vecRikij, Rijr, Rikr, cost, lambda, &
      !               zeta, F1, dF1j, dF1k)
      sf_F1_result = sf_F1_ijk_func(vecRijik, vecRikij, Rijr, Rikr, cost, lambda, zeta)
      
      F1 = sf_F1_result%Fijk
      dF1j = sf_F1_result%dFijk_dRj
      dF1k = sf_F1_result%dFijk_dRk

      Rcr = 1.d0/Rc
      !call sf_F2_ij(Rij, Rc, Rcr, eta, F2ij, dF2ij)
      !call sf_F2_ij(Rik, Rc, Rcr, eta, F2ik, dF2ik)
      !call sf_F2_ij(Rjk, Rc, Rcr, eta, F2jk, dF2jk)

      sf_F2_result = sf_F2_ij_func(Rij, Rc, Rcr, eta)
      F2ij = sf_F2_result(1)
      dF2ij = sf_F2_result(2)

      sf_F2_result = sf_F2_ij_func(Rik, Rc, Rcr, eta)
      F2ik = sf_F2_result(1)
      dF2ik = sf_F2_result(2)

      sf_F2_result = sf_F2_ij_func(Rjk, Rc, Rcr, eta)
      F2jk = sf_F2_result(1)
      dF2jk = sf_F2_result(2)

      G0    = 2.d0*F2ij*F2ik*F2jk
      G4    = F1*G0
      ! factor of 2 for k<j in sf_fingerprint()
!       G(iG) = G(iG) + 2.0d0*G4
      !G(iG) = G(iG) + G4
      G(iG_tmp) = G(iG_tmp) + G4

      if (present(dGi)) then
         F1   = 2.d0*F1
         dG401 = F1*dF2ij*F2ik*F2jk
         dG402 = F1*F2ij*dF2ik*F2jk
         dG403 = F1*F2ij*F2ik*dF2jk
         vecRij0(1:3) = vecRij(1:3)*dG401
         vecRik0(1:3) = vecRik(1:3)*dG402
         vecRjk0(1:3) = vecRjk(1:3)*dG403
         dF1j(1:3) = dF1j(1:3)*G0
         dF1k(1:3) = dF1k(1:3)*G0
         dG4(1:3) = -(dF1j(1:3) + dF1k(1:3)) - vecRij0(1:3) - vecRik0(1:3)
!          dGi(1:3,iG) = dGi(1:3,iG) + 2.0d0*dG4(1:3)
         dGi(1:3,iG_tmp) = dGi(1:3,iG_tmp) + dG4(1:3)
!       end if
!       if (present(dGj)) then
         dG4(1:3) =  dF1j(1:3) + vecRij0(1:3) - vecRjk0(1:3)
!          dGj(1:3,iG) = dGj(1:3,iG) + 2.0d0*dG4(1:3)
         dGj(1:3,iG_tmp) = dGj(1:3,iG_tmp) + dG4(1:3)
!       end if
!       if (present(dGk)) then
         dG4(1:3) =  dF1k(1:3) + vecRik0(1:3) + vecRjk0(1:3)
!          dGk(1:3,iG) = dGk(1:3,iG) + 2.0d0*dG4(1:3)
         dGk(1:3,iG_tmp) = dGk(1:3,iG_tmp) + dG4(1:3)
!       end if
         if (present(dGh)) then
            dGhij(1:3) = ( vecRij0(1:3) + dF1j(1:3) )*Rij
            dGhik(1:3) = ( vecRik0(1:3) + dF1k(1:3) )*Rik
            dGhjk(1:3) =   vecRjk0(1:3)*Rjk
            dGhh(1) = dGhij(1)*vecRij(1)  &
   &                  + dGhik(1)*vecRik(1)  &
   &                  + dGhjk(1)*vecRjk(1)
            dGhh(2) = dGhij(2)*vecRij(2)  &
   &                  + dGhik(2)*vecRik(2)  &
   &                  + dGhjk(2)*vecRjk(2)
            dGhh(3) = dGhij(3)*vecRij(3)  &
   &                  + dGhik(3)*vecRik(3)  &
   &                  + dGhjk(3)*vecRjk(3)
            dGhh(4) = dGhij(2)*vecRij(3)  &
   &                  + dGhik(2)*vecRik(3)  &
   &                  + dGhjk(2)*vecRjk(3)
            dGhh(5) = dGhij(3)*vecRij(1)  &
   &                  + dGhik(3)*vecRik(1)  &
   &                  + dGhjk(3)*vecRjk(1)
            dGhh(6) = dGhij(1)*vecRij(2)  &
   &                  + dGhik(1)*vecRik(2)  &
   &                  + dGhjk(1)*vecRjk(2)
   !           dGh(1:6,iG) = dGh(1:6,iG) + 2.0d0*dGhh(1:6)
            dGh(1:6,iG_tmp) = dGh(1:6,iG_tmp) + dGhh(1:6)
         end if
      end if

    end do
   !omp end target data

    iG = iG + iter

  end subroutine sf_G4_update

  !--------------------------------------------------------------------!
  !         second angular symm. function (eq. 9 in Ref. [1])          !
  !--------------------------------------------------------------------!

  subroutine sf_G5_update(tijk, vecRij, vecRik, Rij, Rik, cost, n, &
                          G, iG, dGi, dGj, dGk, dGh)

    implicit none

    integer,                                    intent(in)    :: tijk
    double precision, dimension(3),             intent(in)    :: vecRij, vecRik
    double precision,                           intent(in)    :: Rij, Rik
    double precision,                           intent(in)    :: cost
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(3,n), optional, intent(inout) :: dGk
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer                        :: iG5
    double precision               :: Rc, Rcr, lambda, zeta, eta
    double precision               :: G5
    double precision, dimension(3) :: dG5, dF1j, dF1k

    double precision :: F1, F2ij, dF2ij, F2ik, dF2ik
    double precision :: Rijr, Rikr, vecRijik(3), vecRikij(3)
    double precision :: dGhh(6), dGhij, dGhik, dF1ij(3), dF1ik(3)

    Rijr = 1.d0/Rij
    Rikr = 1.d0/Rik
    vecRijik(1:3) = -cost*vecRij(1:3) + vecRik(1:3)
    vecRikij(1:3) = -cost*vecRik(1:3) + vecRij(1:3)
    do iG5 = 1, sf_nG3(5,tijk)
       iG     = iG + 1
       Rc     = sf_pG5(1,iG5,tijk)
       lambda = sf_pG5(2,iG5,tijk)
       zeta   = sf_pG5(3,iG5,tijk)
       eta    = sf_pG5(4,iG5,tijk)

       call sf_F1_ijk(vecRijik, vecRikij, Rijr, Rikr, cost, lambda, &
                      zeta, F1, dF1j, dF1k)
       Rcr = 1.d0/Rc
       call sf_F2_ij(Rij, Rc, Rcr, eta, F2ij, dF2ij)
       call sf_F2_ij(Rik, Rc, Rcr, eta, F2ik, dF2ik)

       G5    = F1*F2ij*F2ik
       ! factor of two for k<j in sf_fingerprint()
       G(iG) = G(iG) + 2.0d0*G5

       if (present(dGi)) then
          dG5(1:3) = -(dF1j(1:3)+dF1k(1:3))*F2ij*F2ik        &
                   - vecRij(1:3)*F1*dF2ij*F2ik &
                   - vecRik(1:3)*F1*F2ij*dF2ik
          dGi(1:3,iG) = dGi(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGj)) then
          dG5(1:3) = dF1j(1:3)*F2ij*F2ik        &
                   + vecRij(1:3)*F1*dF2ij*F2ik
          dGj(1:3,iG) = dGj(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGk)) then
          dG5(1:3) = dF1k(1:3)*F2ij*F2ik        &
                   + vecRik(1:3)*F1*F2ij*dF2ik
          dGk(1:3,iG) = dGk(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGh)) then
           dGhij = F1*dF2ij* F2ik*Rij
           dGhik = F1* F2ij*dF2ik*Rik
           dF1ij(1:3) = F2ij*F2ik*dF1j(1:3)*Rij
           dF1ik(1:3) = F2ij*F2ik*dF1k(1:3)*Rik
           dGhh(1) = dGhij*vecRij(1)*vecRij(1)  &
&                  + dGhik*vecRik(1)*vecRik(1)  &
&                  +        dF1ij(1)*vecRij(1)  &
&                  +        dF1ik(1)*vecRik(1)
           dGhh(2) = dGhij*vecRij(2)*vecRij(2)  &
&                  + dGhik*vecRik(2)*vecRik(2)  &
&                  +        dF1ij(2)*vecRij(2)  &
&                  +        dF1ik(2)*vecRik(2)
           dGhh(3) = dGhij*vecRij(3)*vecRij(3)  &
&                  + dGhik*vecRik(3)*vecRik(3)  &
&                  +        dF1ij(3)*vecRij(3)  &
&                  +        dF1ik(3)*vecRik(3)
           dGhh(4) = dGhij*vecRij(2)*vecRij(3)  &
&                  + dGhik*vecRik(2)*vecRik(3)  &
&                  +        dF1ij(2)*vecRij(3)  &
&                  +        dF1ik(2)*vecRik(3)
           dGhh(5) = dGhij*vecRij(3)*vecRij(1)  &
&                  + dGhik*vecRik(3)*vecRik(1)  &
&                  +        dF1ij(3)*vecRij(1)  &
&                  +        dF1ik(3)*vecRik(1)
           dGhh(6) = dGhij*vecRij(1)*vecRij(2)  &
&                  + dGhik*vecRik(1)*vecRik(2)  &
&                  +        dF1ij(1)*vecRij(2)  &
&                  +        dF1ik(1)*vecRik(2)
           dGh(1:6,iG) = dGh(1:6,iG) + 2.0d0*dGhh(1:6)
       end if

    end do

  end subroutine sf_G5_update


  !======================= auxiliary functions ========================!



  !--------------------------------------------------------------------!
  !                      cosine cut-off function                       !
  !                        (eq. 4 in Ref. [1])                         !
  !--------------------------------------------------------------------!

  subroutine sf_cut(Rij, Rc, Rcr, fc, dfc)

    implicit none

    double precision, intent(in)  :: Rij, Rc, Rcr
    double precision, intent(out) :: fc, dfc
    real(8) :: x, d
    integer :: m

    if (Rij >= Rc) then
       fc  = 0.0d0
       dfc = 0.0d0
       return
    end if

    fc  =  0.5d0*(cos(PI*Rij/Rc) + 1.0d0)
    dfc = -0.5d0*PI/Rc*sin(PI*Rij/Rc)

    return  ! kn


    x = PI*Rij*Rcr
    if( norder == 1 ) then
        !---first-order interpolation
        d = x*dxr
        m = d
        d = d - m
         fc =  (1d0-d)*cutv(m)+d*cutv(m+1)
        dfc = ((1d0-d)*cutd(m)+d*cutd(m+1))*Rcr
    else
        !---second-order interpolation
        m = x/(2.d0*dx)
        m = 2*m
        d = 0.5d0*( x*dxr - dble(m) )
         fc =  d*( (d-1.d0)*cutva(m) + cutv(m+2) - cutv(m) ) + cutv(m)
        dfc = (d*( (d-1.d0)*cutda(m) + cutd(m+2) - cutd(m) ) + cutd(m))*Rcr
    end if

  end subroutine sf_cut


  subroutine sf_set_table

    implicit none
    integer :: k
    real(8) :: x

    if( norder == 1 ) then
        nbin = 1000000
    else
        nbin = 10000
    end if

    allocate( cutv(0:nbin), cutd(0:nbin), exv(0:nbin) )

!    pi = acos(-1.d0)

    !---cos, sin for sf_cut
    dx  = pi/(nbin-1)
    dxr = 1.d0/dx

     cutv(:) = 0.d0
     cutd(:) = 0.d0
    do k = 0, nbin
       x = dble(k)*dx
       cutv(k) = 0.5d0*(cos(x) + 1.d0)
       cutd(k) = -0.5d0*pi*sin(x)
    end do

    !---exp
    exx = 40.d0
    dex  = exx/(nbin-1)
    dexr = 1.d0/dex
    
    exv(:) = 0.d0
    do k = 0, nbin
       x = dble(k)*dex
       exv(k) = exp(-x)
    end do


    if( norder /= 1 ) then
        !---for the second-order interpolation
        allocate( cutva(0:nbin), cutda(0:nbin), exva(0:nbin) )

        cutva(:) = 0.d0
        cutda(:) = 0.d0
         exva(:) = 0.d0
        do k = 0, nbin - 2
           cutva(k) = 2.d0*( cutv(k) + cutv(k+2) - 2.d0*cutv(k+1) )
           cutda(k) = 2.d0*( cutd(k) + cutd(k+2) - 2.d0*cutd(k+1) )
            exva(k) = 2.d0*(  exv(k) +  exv(k+2) - 2.d0* exv(k+1) )
        end do
    end if


  end subroutine sf_set_table

end module symmfunc

