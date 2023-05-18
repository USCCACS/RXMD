module stat_mod

  implicit none
  real(8),parameter :: pi = 4.d0*datan(1.d0)
  integer,parameter :: NTABLES = 200, NTABLES_BA=180
  real(8),parameter :: RCUT = 6d0, DRI = NTABLES/RCUT
  !real(8),parameter :: RCUT = 15d0, DRI = NTABLES/RCUT
  real(8),parameter :: QCUT = 10d0, DQ = QCUT/NTABLES
  integer,parameter :: MAXNEIGHBS = 500

  type bond_length 
     character(len=2) :: A,B
     real(8) :: rc
  end type

  real(8),parameter :: BARC0=4d0
  type(bond_length) :: bond_length0(2) = [ &
                       bond_length('H','N',1.2d0), &
                       bond_length('Na','O',3.0d0) &
                       ]
                       !bond_length('H','O',1.0d0), &
                       !bond_length('Na','O',2.0d0) &

  type NSD_type ! Neutron Scattering Data type
     character(len=2) :: name
     real(8) :: length
  end type

  ! neutron scattering length data are from 
  !  https://www.nist.gov/ncnr/neutron-scattering-lengths-list
  type(NSD_type),parameter :: NSD0(14)=[&
          NSD_type(name='Ge',length=8.185d-5),  NSD_type(name='Se',length=7.970d-5), &
          NSD_type(name='Sb',length=5.57d-5),   NSD_type(name='Te',length=5.80d-5), & 
          NSD_type(name='C', length=6.646d-5),  NSD_type(name='Si',length=4.1491d-5), &
          NSD_type(name='O', length=5.803d-5),  NSD_type(name='Al',length=3.449d-5), &
          NSD_type(name='H', length=-3.7390d-5),NSD_type(name='Na',length=3.63d-5), &
          NSD_type(name='Cl',length=9.5770d-5), NSD_type(name='Pb',length=9.405d-5), &
          NSD_type(name='Ti',length=-3.4380d-5),NSD_type(name='N', length=9.36d-5)  &
          ]

  type AFF_type ! Atomic Form Factor type
    character(len=2) :: name
    real(8) :: c
    real(8) :: a(4)
    real(8) :: b(4)
  end type

  ! atomic form factor data are taken from
  ! http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
  type(AFF_type),parameter :: AFF0(7)=[ &
          AFF_type(name='H', a=[0.489918d0,0.262003d0,0.196767d0,0.049879d0], &
                   b=[20.6593d0,7.74039d0,49.5519d0,2.20159d0], c=0.001305d0 ), &
          AFF_type(name='O', a=[3.0485d0, 2.2868d0, 1.5463d0, 0.867d0], &
                   b=[13.2771d0, 5.7011d0, 0.3239d0, 32.9089d0], c=0.2508d0), &
          AFF_type(name='C', a=[2.31d0, 1.02d0, 1.5886d0, 0.865d0], &
                   b=[20.8439d0, 10.2075d0, 0.5687d0, 51.6512d0], c=0.2156d0), &
          AFF_type(name='N', a=[12.2126d0, 3.1322d0, 2.0125d0, 1.1663d0], &
                   b=[0.0057d0, 9.8933d0, 28.9975d0, 0.5826d0], c=-11.529d0), &
          AFF_type(name='Li',a=[1.1282d0, 0.7508d0, 0.6175d0, 0.4653d0], &
                   b=[3.9546d0, 1.0524d0, 85.3905d0, 168.261d0], c=0.0377d0), &
          AFF_type(name='Na',a=[4.7626d0, 3.1736d0, 1.2674d0, 1.1128d0], &
                   b=[3.285d0, 8.8422d0, 0.3136d0, 129.424d0], c=0.676), &
          AFF_type(name='K', a=[8.2186d0, 7.4398d0, 1.0519d0, 0.8659d0], &
                   b=[12.7949d0, 0.7748d0, 213.187d0, 41.6841d0], c=1.4228d0) &
          ]


  type base_atom_type
    real(8) :: pos(3), rr
    integer :: itype, ir, id
    character(len=:),allocatable :: elem
  end type

  type, extends(base_atom_type) :: nbrlist_type
    type(base_atom_type),allocatable :: nbrs(:)
    integer :: counter
  end type

  type string_array
     character(len=:),allocatable :: str
  end type

  type mdframe
     character(len=:),allocatable :: filename
     type(string_array),allocatable :: elems(:)
     real(8),allocatable :: pos(:,:), v(:,:), f(:,:), q(:)
     integer,allocatable :: itype(:)
     real(8) :: lattice(6)
     integer :: num_atoms
  end type

  type analysis_context

     type(string_array),allocatable :: elems(:)
     type(NSD_type),allocatable :: NSD(:)
     type(AFF_type),allocatable :: AFF(:)

     real(8),allocatable :: concentration(:), gr(:,:,:), nr(:,:,:), ba(:,:,:,:), sq(:,:,:)

     real(8),allocatable :: ba_rc(:,:)

     integer :: num_atoms, num_atom_types
     real(8) :: volume, lattice(6), kvector(3)
     integer,allocatable :: num_atoms_per_type(:)

     integer :: num_sample_frames 

  contains

     procedure :: print => print_analysis_context 
     procedure :: save_stat => save_analysis_results

  end type

contains

!-----------------------------------------------------------------------------------------
  function xform(aff, qval) result(fx) 
!-----------------------------------------------------------------------------------------
      type(AFF_type) :: aff
      real(8) :: pi4_i = 1d0/(4d0*pi)
      real(8) :: fx, qval
      integer :: i
      fx = aff%c
      do i = 1, size(aff%a)
         fx = fx + aff%a(i)*exp(-aff%b(i)*(qval*pi4_i)**2)
      enddo
  end function

!-----------------------------------------------------------------------------------------
  function nbrlist_ctor_from_mdframe(oneframe) result(c)
!-----------------------------------------------------------------------------------------
     type(mdframe),intent(in) :: oneframe
     type(nbrlist_type),allocatable :: c(:)
     integer :: i

     allocate(c(oneframe%num_atoms))

     do i=1, oneframe%num_atoms
       c(i)%counter = 0

       c(i)%pos(1:3) = oneframe%pos(i,1:3)
       c(i)%itype = oneframe%itype(i)
       c(i)%elem = oneframe%elems(i)%str
       c(i)%id = i
!print*,c(i)%pos(1:3), c(i)%itype, c(i)%elem , c(i)%id

       allocate(c(i)%nbrs(MAXNEIGHBS)) !FIXME automatically find the size
     enddo

!print*,repeat('-=',60)
     return
  end function

!-----------------------------------------------------------------------------------------
  subroutine print_analysis_context(this) 
!-----------------------------------------------------------------------------------------
     class(analysis_context) :: this
     character(len=:),allocatable :: name
     integer :: i,j

     print'(a)',repeat('-',60)
     print'(a,i6)','num_atoms: ', this%num_atoms
     print'(a,3f8.3, 3f8.2)','lattice: ', this%lattice
     print'(a,es15.5)','volume: ', this%volume
     print'(a,3f8.3)','k-vector: ', this%kvector

     print'(a, i6)', 'num_elements: ', size(this%elems)
     print'(a $)', 'elements: '
     do i = 1, size(this%elems)
        name = this%elems(i)%str
        write(6, fmt='(i3, a)', advance='no') get_index(this%elems, name),'-'//adjustl(name)//'   '
     enddo
     print*
     print'(a,10i6)', 'num_atoms_per_type: ', this%num_atoms_per_type
     print'(a,10f8.5)', 'concentration: ', this%concentration

     print'(a)',repeat('-',60)
     print'(a $)', 'neutron scattering length: '
     do i=1, size(this%NSD)
         write(6,fmt='(i3,a3,es10.3,a2)', advance='no')  i, '-'//this%NSD(i)%name, this%NSD(i)%length, ', '
     enddo
     print*

     print'(a)', 'x-ray form factor: '
     do i=1, size(this%AFF)
        print'(a $)','a1,b1,a2,b2,a3,b3,a4,b4,c '//this%AFF(i)%name
        do j=1, size(this%AFF(i)%a)
           print'(2f12.5 $)', this%AFF(i)%a(j),this%AFF(i)%b(j)
        enddo
        print'(f12.5)',this%AFF(i)%c
     enddo

     print'(a)',repeat('-',60)

  end subroutine

!-----------------------------------------------------------------------------------------
  subroutine save_analysis_results(this) 
!-----------------------------------------------------------------------------------------
     class(analysis_context) :: this
     integer :: iunit,ity,jty,kty,k,kk,l
     real(8) :: dr, rho, dqk, Sxq, Sxq_denom, Snq, Snq_denom, Gnr, Gnr_denom, prefactor, prefactor2, bavalue
     character(len=1) :: a1

     ! get the number density, rho
     rho = this%num_atoms/this%volume

     open(newunit=iunit,file='gr.dat',form='formatted')

     ! g(r) header part
     write(unit=iunit,fmt='(a $)') ' 1-distance,'
     do ity=1,size(this%elems); 
     do jty=1,size(this%elems)
        write(unit=iunit,fmt='(a12, 1x)', advance='no') & 
           this%elems(ity)%str//'-'//this%elems(jty)%str//'(gr),'
     enddo; enddo
     do ity=1,size(this%elems);
     do jty=1,size(this%elems)
        write(unit=iunit,fmt='(a12, 1x)', advance='no') & 
           this%elems(ity)%str//'-'//this%elems(jty)%str//'(nr),'
     enddo; enddo
     write(unit=iunit,fmt='(a)') ' neutron_gr'

     Gnr_denom = sum(this%NSD(:)%length * this%concentration(:))
     Gnr_denom = Gnr_denom**2

     ! get coordination number, n(r)
     do k=1,size(this%gr,dim=3)
        do ity=1,size(this%elems)
        do jty=1,size(this%elems)
          this%nr(ity,jty,k) = sum(this%gr(ity,jty,1:k))/(this%num_atoms_per_type(ity)*this%num_sample_frames)
        enddo; enddo
     enddo

     do k=1,size(this%gr,dim=3)

        dr = k/DRI
        prefactor = 4.d0*pi*dr*dr*rho/DRI
        write(unit=iunit,fmt='(f12.5,1x,a1)',advance='no') dr,','

        Gnr = 0.d0
        do ity=1,size(this%elems)
        do jty=1,size(this%elems)

          prefactor2 = this%concentration(jty) * this%num_atoms_per_type(ity) * this%num_sample_frames
          this%gr(ity,jty,k) = this%gr(ity,jty,k)/(prefactor*prefactor2)
          write(iunit, fmt='(f12.5,1x,1a)',advance='no') this%gr(ity,jty,k),','

          Gnr = Gnr + this%gr(ity,jty,k) * &
                      this%concentration(ity) * this%concentration(jty) * & 
                      this%NSD(ity)%length * this%NSD(jty)%length
        enddo; enddo

        do ity=1,size(this%elems)
        do jty=1,size(this%elems)
           write(iunit, fmt='(f12.5,1x,a1)',advance='no') this%nr(ity,jty,k),','
        enddo; enddo
        write(iunit, fmt='(f12.5)') Gnr/Gnr_denom
     enddo

     close(iunit)

     open(newunit=iunit,file='sq.dat',form='formatted')
     write(unit=iunit,fmt='(a)',advance='no') ' wave_number,'
     do ity=1,size(this%elems)
     do jty=1,size(this%elems)
        write(unit=iunit,fmt='(a5,a5,a1)',advance='no') & 
           '     ',this%elems(ity)%str//'-'//this%elems(jty)%str,','
     enddo; enddo
     write(unit=iunit,fmt='(a)') '  Snq,   Sxq'

     ! get the denominator of Sn(q)
     Snq_denom = sum(this%NSD(:)%length * this%concentration(:))
     Snq_denom = Snq_denom**2

     do kk=1, size(this%sq,dim=3)

        dqk = DQ*kk

        ! the denominator for Sx(q) should be q-dependent.
        Sxq_denom = 0.d0
        do ity=1, size(this%elems)
           Sxq_denom = Sxq_denom + xform(this%AFF(ity),dqk) * this%concentration(ity)
        enddo 
        Sxq_denom = Sxq_denom**2

        write(unit=iunit,fmt='(f12.5,a1)',advance='no') dqk,','

        Snq = 0.d0; Sxq = 0.d0
        do ity=1,size(this%elems)
        do jty=1,size(this%elems)

           ! initialize sq for the r-value integral 
           this%sq(ity,jty,kk) = 0.d0

           do k=1,size(this%gr,dim=3)
              dr = k/DRI
              !prefactor = (this%gr(ity,jty,k)-1.d0)*cos(0.5d0*pi*dr/RCUT) ! (gr - 1)*taper
              prefactor = (this%gr(ity,jty,k)-1.d0) ! (gr - 1)*taper
              this%sq(ity,jty,kk) = this%sq(ity,jty,kk) + &
                  dr*dr*prefactor*sin(dr*dqk)/(dr*dqk)/DRI
           enddo

           ! get prefactor 
           prefactor = 4d0*pi*rho*sqrt(this%concentration(ity)*this%concentration(jty))
           this%sq(ity,jty,kk) = this%sq(ity,jty,kk)*prefactor

           ! add delta function if ity==jty
           if(ity==jty) this%sq(ity,jty,kk) = this%sq(ity,jty,kk) + 1.d0
           write(iunit, fmt='(f12.5,1x,a1)',advance='no') this%sq(ity,jty,kk), ','

           Snq = Snq + this%sq(ity,jty,kk) * &
                       sqrt(this%concentration(ity) * this%concentration(jty)) * & 
                       this%NSD(ity)%length * this%NSD(jty)%length

           Sxq = Sxq + this%sq(ity,jty,kk) * &
                       sqrt(this%concentration(ity) * this%concentration(jty)) * & 
                       xform(this%AFF(ity),dqk) * xform(this%AFF(jty),dqk)
        enddo; enddo
        write(iunit, fmt='(f12.5,a2,f12.5)') Snq/Snq_denom, ',', Sxq/Sxq_denom

     enddo
     close(iunit)

     open(newunit=iunit,file='ba.dat',form='formatted')

     write(unit=iunit,fmt='(a,1x)',advance='no') '1-angle,'
     do ity=1,size(this%elems)
     do jty=1,size(this%elems)
     do kty=1,size(this%elems)
        write(unit=iunit,fmt='(a12,a1,1x)',advance='no') & 
           this%elems(jty)%str//'-'//this%elems(ity)%str//'-'//this%elems(kty)%str,','
     enddo; enddo; enddo
     write(unit=iunit,fmt=*)

     do k=1,size(this%ba,dim=4)
        write(unit=iunit,fmt='(i8,a1)',advance='no') k,','
        do ity=1,size(this%elems)
        do jty=1,size(this%elems)
        do kty=1,size(this%elems)
           bavalue = sum(this%ba(jty,ity,kty,:))
           if(bavalue>0.d0) then
             write(unit=iunit,fmt='(es12.5,a1,1x)',advance='no') &
               this%ba(jty,ity,kty,k)/(bavalue*this%num_atoms_per_type(ity)*this%num_sample_frames),','
           else
             write(unit=iunit,fmt='(es12.5,a1,1x)',advance='no') 0.d0,','
           endif
        enddo; enddo; enddo
        write(unit=iunit,fmt=*)
     enddo

     close(iunit)

  end subroutine


!-----------------------------------------------------------------------------------------
  function get_analysis_context_from_mdframe(oneframe) result(c)
!-----------------------------------------------------------------------------------------
     type(mdframe),intent(in) :: oneframe
     type(analysis_context) :: c
     integer :: i, j, ne, ity, jty, kty, idx

     character(len=256) :: argv
     integer :: iunit

     character(len=2) :: A,B

     c%num_atoms = oneframe%num_atoms
     c%lattice = oneframe%lattice
     c%volume = oneframe%lattice(1)*oneframe%lattice(2)*oneframe%lattice(3)

     allocate(c%elems(0))
     do i=1, size(oneframe%elems)
        if( get_index(c%elems,oneframe%elems(i)%str) < 0 ) &
            c%elems = [c%elems, string_array(oneframe%elems(i)%str)]
     enddo

     ! setup neutron scattering length data for existing atom types
     allocate(c%NSD(0), c%AFF(0))
     do i=1, size(c%elems)
        do j=1, size(NSD0)
           if(c%elems(i)%str==NSD0(j)%name) c%NSD = [c%NSD, NSD0(j)]
        enddo
        do j=1, size(AFF0)
           if(c%elems(i)%str==AFF0(j)%name) c%AFF = [c%AFF, AFF0(j)]
        enddo
     enddo

     ! check if NSD is found for all elements
     if ( size(c%NSD) /= size(c%elems) ) then
        print*,'Error : missing NSD', size(c%NSD), c%NSD
        stop
     endif

     ! check if AFF is found for all elements
     if ( size(c%AFF) /= size(c%elems) ) then
        print*,'Error : missing AFF', size(c%AFF), c%AFF
        stop
     endif

     c%kvector(1:3)=2d0*pi/c%lattice(1:3)

     ne = size(c%elems)
     c%num_atom_types = ne

     allocate(c%gr(ne,ne,NTABLES),c%nr(ne,ne,NTABLES),c%sq(ne,ne,NTABLES),c%ba(ne,ne,ne,NTABLES_BA))
     c%gr=0.d0;  c%nr=0.d0;  c%sq=0.d0; c%ba=0.d0

     allocate(c%concentration(ne), c%num_atoms_per_type(ne))
     c%concentration=0.0d0; c%num_atoms_per_type=0

     do i=1, size(oneframe%elems)
        idx = get_index(c%elems, oneframe%elems(i)%str) 
        if(idx>0) then 
            c%num_atoms_per_type(idx) = c%num_atoms_per_type(idx)+1
        else
            print'(a,a,i6)', 'ERROR: uncategorized atom found: in ', oneframe%filename, i
        endif
     enddo

     c%concentration = dble(c%num_atoms_per_type)/c%num_atoms

     c%num_sample_frames = 0

     print'(a)',repeat('-',60)
     ! setup bond angle cutoff
     allocate(c%ba_rc(ne,ne))
     c%ba_rc=BARC0
     do ity = 1, ne
     do jty = 1, ne
        do i=1,size(bond_length0)
           if( bond_length0(i)%A ==  c%elems(jty)%str .and.  bond_length0(i)%B == c%elems(ity)%str) then
              c%ba_rc(ity,jty) = bond_length0(i)%rc
              c%ba_rc(jty,ity) = bond_length0(i)%rc
              print'(a,f6.2)',' found bond angle cutff '//trim(bond_length0(i)%A)//'-'//trim(bond_length0(i)%B), bond_length0(i)%rc
           endif
        enddo
     enddo
     enddo
     print'(a)',repeat('-',60)

     call c%print()
  end function

!-----------------------------------------------------------------------------------------
  function get_index(elems, name) result(idx)
!-----------------------------------------------------------------------------------------
    type(string_array),allocatable,intent(in) :: elems(:)
    character(len=:),allocatable,intent(in) :: name
    integer :: idx

    do idx = 1, size(elems)
       if (elems(idx)%str==name) return 
    enddo
    idx = -1

    return
  end

end module
