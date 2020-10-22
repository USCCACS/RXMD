module stat_mod

  implicit none
  real(8),parameter :: pi = 4.d0*datan(1.d0)
  integer,parameter :: NTABLES = 50000, NTABLES_SQ=1000, NTABLES_BA=180
  real(8),parameter :: RCUT = 3d0, DRI = NTABLES/RCUT
  !real(8),parameter :: RCUT = 15d0, DRI = NTABLES/RCUT
  real(8),parameter :: QCUT = 6d0, DQ = QCUT/NTABLES_SQ
  integer,parameter :: MAXNEIGHBS = 100

  type NSD_type ! Neutron Scattering Data type
     character(len=2) :: name
     real(8) :: length
  end type

  ! neutron scattering length data are from 
  !  https://www.nist.gov/ncnr/neutron-scattering-lengths-list
  type(NSD_type) :: NSD0(11)=[NSD_type(name='Ge',length=8.185d-5), & 
                              NSD_type(name='Se',length=7.970d-5), &
                              NSD_type(name='Sb',length=5.57d-5), &
                              NSD_type(name='Te',length=5.80d-5), & 
                              NSD_type(name='C',length=6.646d-5), &
                              NSD_type(name='Si',length=4.1491d-5), &
                              NSD_type(name='O',length=5.803d-5), &
                              NSD_type(name='Al',length=3.449d-5), &
                              NSD_type(name='H',length=-3.7390d-5), &
                              NSD_type(name='Na',length=3.63d-5), &
                              NSD_type(name='Cl',length=9.5770d-5) ]

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

  integer,parameter :: NUM_BA=1
  real(8),parameter :: BA_CUTOFF=2.0d0

  type analysis_context

     type(string_array),allocatable :: elems(:)
     type(NSD_type),allocatable :: NSD(:)

     real(8),allocatable :: concentration(:), gr(:,:,:), nr(:,:,:), ba(:,:,:,:,:), sq(:,:,:)
     real(8) :: rc_ba(NUM_BA)

     integer :: num_atoms, num_atom_types
     real(8) :: volume, lattice(6), kvector(3)
     integer,allocatable :: num_atoms_per_type(:)

     real(8) :: box_rc, box_dri

     integer :: num_sample_frames 

  contains

     procedure :: print => print_analysis_context 
     procedure :: save_stat => save_analysis_results

  end type

contains

!-----------------------------------------------------------------------------------------
  function base_atom_type_ctor(pos_, itype_, rr_, ir_, id_, elem_) result(a)
!-----------------------------------------------------------------------------------------
  type(base_atom_type) :: a
    real(8),optional :: pos_(3), rr_
    integer,optional :: itype_, ir_, id_
    character(len=:),optional,allocatable :: elem_

    a%pos = 0.d0; a%rr = 0.d0
    a%itype = 0;  a%ir = 0;  a%id = 0
    a%elem = "NA"

    if(present(pos_)) a%pos = pos_ 
    if(present(itype_)) a%itype = itype_
    if(present(rr_)) a%rr = rr_
    if(present(ir_)) a%ir = ir_
    if(present(id_)) a%id = id_
    if(present(elem_)) a%elem = elem_

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
     integer :: i

     print'(a)',repeat('-',60)
     print'(a,i12)','num_atoms: ', this%num_atoms
     print'(a,3f8.3, 3f8.2)','lattice: ', this%lattice
     print'(a,es15.5)','volume: ', this%volume
     print'(a,3f8.3)','k-vector: ', this%kvector

     print'(a, i6)', 'num_elements: ', size(this%elems)
     write(6,fmt='(a)', advance='no') 'elements: '
     do i = 1, size(this%elems)
        name = this%elems(i)%str
        write(6, fmt='(i3, a)', advance='no') get_index(this%elems, name),'-'//adjustl(name)//'   '
     enddo
     print*
     print'(a,10i9)', 'num_atoms_per_type: ', this%num_atoms_per_type
     print'(a,10f8.5)', 'concentration: ', this%concentration
     write(6,fmt='(a)',advance='no') 'neutron scattering: '
     do i=1, size(this%NSD)
         write(6,fmt='(i3,a3,es10.3,a2)', advance='no')  i, '-'//this%NSD(i)%name, this%NSD(i)%length, ', '
     enddo
     print*
     print'(a,2f8.3)', 'box_rc,box_dri: ', this%box_rc, this%box_dri
     print'(a)',repeat('-',60)

  end subroutine

!-----------------------------------------------------------------------------------------
  subroutine save_analysis_results(this) 
!-----------------------------------------------------------------------------------------
     class(analysis_context) :: this

     integer :: iunit,ity,jty,kty,k,kk,l
     real(8) :: dr, rho, dqk, Snq, Snq_denom, Gnr, Gnr_denom, prefactor, prefactor2, bavalue
     character(len=1) :: a1

     ! get the number density, rho
     rho = this%num_atoms/this%volume

     open(newunit=iunit,file='gr.dat',form='formatted')

     ! g(r) header part
     write(unit=iunit,fmt='(a)',advance='no') ' 1-distance,'
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

        dr = k/this%box_dri
        prefactor = 4.d0*pi*dr*dr*rho/this%box_dri
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
     write(unit=iunit,fmt='(a)') '  Snq'

     ! get the denominator of Sn(q) 
     !Snq_denom = sum(this%NSD(:)%length * this%concentration(:))
     Snq_denom = sum(this%concentration(:))
     Snq_denom = Snq_denom**2

     do kk=1, size(this%sq,dim=3)

        dqk = DQ*kk
        write(unit=iunit,fmt='(f12.5,a1)',advance='no') dqk,','

        Snq = 0.d0
        do ity=1,size(this%elems)
        do jty=1,size(this%elems)

           ! initialize sq for the r-value integral 
           this%sq(ity,jty,kk) = 0.d0

           do k=1,size(this%gr,dim=3)
              dr = k/this%box_dri
              prefactor = (this%gr(ity,jty,k)-1.d0)*cos(0.5d0*pi*dr/this%box_rc) ! (gr - 1)*taper
              !prefactor = (this%gr(ity,jty,k)-1.d0) ! (gr - 1)*taper
              this%sq(ity,jty,kk) = this%sq(ity,jty,kk) + &
                  dr*dr*prefactor*sin(dr*dqk)/(dr*dqk)/this%box_dri
           enddo

           ! get prefactor 
           prefactor = 4d0*pi*rho*sqrt(this%concentration(ity)*this%concentration(jty))
           this%sq(ity,jty,kk) = this%sq(ity,jty,kk)*prefactor

           ! add delta function if ity==jty
           if(ity==jty) this%sq(ity,jty,kk) = this%sq(ity,jty,kk) + 1.d0
           write(iunit, fmt='(f12.5,1x,a1)',advance='no') this%sq(ity,jty,kk), ','

           !Snq = Snq + this%sq(ity,jty,kk) * &
           !            sqrt(this%concentration(ity) * this%concentration(jty)) * & 
           !            this%NSD(ity)%length * this%NSD(jty)%length
           Snq = Snq + this%sq(ity,jty,kk) * &
                       sqrt(this%concentration(ity) * this%concentration(jty)) 
        enddo; enddo
        write(iunit, fmt='(f12.5)') Snq/Snq_denom

     enddo
     close(iunit)

     do ity=1,size(this%elems)
        open(newunit=iunit,file='ba-'//this%elems(ity)%str//'.dat',form='formatted')

        write(unit=iunit,fmt='(a,1x)',advance='no') '1-angle,'
        do jty=1,size(this%elems)
        do kty=1,size(this%elems)
           do l=1,size(this%ba,dim=4)
             write(a1,'(i1)') l
               write(unit=iunit,fmt='(a12,a1,1x)',advance='no') & 
                 this%elems(jty)%str//'-'//this%elems(ity)%str//'-'//this%elems(kty)%str//'_'//a1,','
           enddo
        enddo; enddo
        write(unit=iunit,fmt=*)

        do k=1,size(this%ba,dim=5)
           write(unit=iunit,fmt='(i8,a1)',advance='no') k,','
           do jty=1,size(this%elems)
           do kty=1,size(this%elems)
              do l=1,size(this%ba,dim=4)
                 bavalue = sum(this%ba(jty,ity,kty,l,:))
                 if(bavalue>0.d0) then
                   write(unit=iunit,fmt='(es12.5,a1,1x)',advance='no') &
                     this%ba(jty,ity,kty,l,k)/(bavalue*this%num_atoms_per_type(ity)*this%num_sample_frames),','
                 else
                   write(unit=iunit,fmt='(es12.5,a1,1x)',advance='no') 0.d0,','
                 endif
              enddo
           enddo; enddo
           write(unit=iunit,fmt=*)
        enddo

        close(iunit)
     enddo

  end subroutine


!-----------------------------------------------------------------------------------------
  function get_analysis_context_from_mdframe(oneframe) result(c)
!-----------------------------------------------------------------------------------------
     type(mdframe),intent(in) :: oneframe
     type(analysis_context) :: c
     integer :: i, j, ne, ity, jty, kty, idx, num_elems0
     type(string_array),allocatable :: elems0(:)

     c%num_atoms = oneframe%num_atoms
     c%lattice = oneframe%lattice
     c%volume = oneframe%lattice(1)*oneframe%lattice(2)*oneframe%lattice(3)


     allocate(elems0(10)); num_elems0=0
     do i=1, size(elems0); elems0(i)%str = ''; enddo

     do i=1, size(oneframe%elems)
        !if( get_index(c%elems,oneframe%elems(i)%str) < 0 ) &
        !    c%elems = [c%elems, string_array(oneframe%elems(i)%str)]
        if( get_index(elems0,oneframe%elems(i)%str) < 0 ) then
           num_elems0 = num_elems0 + 1
           elems0(num_elems0) = string_array(oneframe%elems(i)%str)
        endif
     enddo
     allocate(c%elems(num_elems0))
     do i=1, num_elems0
        c%elems(i)%str = elems0(i)%str
     enddo
    

     ! setup neutron scattering length data for existing atom types
     allocate(c%NSD(size(c%elems)))
     do i=1, size(c%elems)
        do j=1, size(NSD0)
           if(c%elems(i)%str==NSD0(j)%name) c%NSD(i) = NSD0(j)
        enddo
     enddo

     c%kvector(1:3)=2d0*pi/c%lattice(1:3)

     ne = size(c%elems)
     c%num_atom_types = ne

     allocate(c%gr(ne,ne,NTABLES),c%nr(ne,ne,NTABLES),c%sq(ne,ne,NTABLES_SQ),c%ba(ne,ne,ne,NUM_BA,NTABLES_BA))
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

     do i=1, NUM_BA
       c%rc_ba(i)=BA_CUTOFF*dble(i)/NUM_BA
     enddo

     c%box_rc = minval(oneframe%lattice(1:3))*0.5d0
     c%box_dri= size(c%gr,dim=3)/c%box_rc

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