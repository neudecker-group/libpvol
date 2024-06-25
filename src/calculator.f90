!================================================================================!
! This file is part of xhcfflib.
!
! Copyright (C) 2023 Felix Zeller, Tim Neudecker, Philipp Pracht
!
! xhcfflib is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xhcfflib is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xhcfflib.  If not, see <https://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------!
!> Adapted from the xtb GBSA source code which can be found under the
!> GNU LGPL 3.0 license, Copyright (C) 2019-2020 Sebastian Ehlert
!> at https://github.com/grimme-lab/xtb
!================================================================================!

!> Implementation of the libpv calculator
!> this is mostly boiler plate code to set up the grid tesselation.
!> the evaluation of the energy and gradient are carried out during the update call
!> grad and energycalcs are implemented in the engrad routine respective to the model
module libpv_surface_engine
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use xhcff_surface_math_wrapper,only:matDet3x3,dot,gemv,symv
  use xhcff_surface_search,only:bisectSearch
  use xhcff_surface_lebedev,only:gridSize,getAngGrid
  use pv_engrad, only:pv_eg
  use xhcff_engrad, only:xhcff_eg
  use xhcff_surface_vdwradd3,only:vanDerWaalsRadD3,vanDerWaalsRadBondi

  implicit none
  private

  public :: surface_engine
  public :: lebedev

  type :: surface_engine

    !> number of atoms
    integer :: nat

    !> atom types
    integer,allocatable :: at(:)

    !> number of pairs
    integer :: ntpair

    !> number of angular grid points
    integer :: nang

    !> angular grid
    real(wp),allocatable :: angGrid(:,:)
    real(wp),allocatable :: angWeight(:)

    !> van der Waals radii of the particles
    real(wp),allocatable :: vdwr(:)

    !> pair descreening approximation radii
    real(wp),allocatable :: rho(:)

    !> offset van der Waals radii
    real(wp),allocatable :: svdw(:)

    !> cut-off radius for the Born radius NN list
    real(wp) :: lrcut

    !> cut-off radius for the SASA NN list
    real(wp) :: srcut

    !> number of neighbors for Born radii
    integer :: nnrad

    !> number of neighbors for SASA computation
    integer,allocatable :: nnsas(:)

    !> neighbors of an atom for Born radii
    integer,allocatable :: nnlistr(:,:)

    !> neighbors of an atom for SASA
    integer,allocatable :: nnlists(:,:)

    !> all pairs indeces array
    integer,allocatable :: ppind(:,:)

    !> all pairs vector differences and magnitudes array
    real(wp),allocatable :: ddpair(:,:)

    !> Atom specific surface data
    real(wp),allocatable :: vdwsa(:)
    real(wp),allocatable :: wrp(:)
    real(wp),allocatable :: trj2(:,:)

    !> Atomic surfaces
    real(wp) :: gsasa
    real(wp) :: sasagam
    real(wp),allocatable :: gamsasa(:)
    real(wp) :: sasa

    !> Molecular pressure gradient
    real(wp), allocatable :: grad(:,:)
    real(wp) :: energy
    real(wp) :: surface
    real(wp) :: volume
    real(wp) :: pressure

    !> pressure model
    !> 0: PV
    !> 1  XHCFF
    integer :: model

    !> Probe radius in Bohr and Angstroem
    real(wp) :: probeRad_aa
    real(wp) :: probeRad_au

  contains

    !> Update coordinates and internal state
    procedure :: update

    !> print
    procedure :: info

    !> initialization
    procedure :: init => init_engine

    !> setup (init+update)
    procedure :: setup => setup_engine

    !> deallocator
    procedure :: deallocate => deallocate_engine

  end type surface_engine

  real(wp),parameter :: autoaa = 0.52917726_wp
  real(wp),parameter :: aatoau = 1.0_wp/autoaa
  real(wp),parameter :: autokcal = 627.50947428_wp
  real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
  real(wp),parameter :: fourpi = 4.0_wp*pi

  !> Smoothing dielectric function parameters
  real(wp),parameter :: w = 0.3_wp*aatoau
  real(wp),parameter :: w3 = w*(w*w)
  real(wp),parameter :: ah0 = 0.5_wp
  real(wp),parameter :: ah1 = 3._wp/(4.0_wp*w)
  real(wp),parameter :: ah3 = -1._wp/(4.0_wp*w3)

  !> Surface calculation cutoffs
  real(wp),parameter :: lrcut = 35.0_wp*aatoau
  real(wp),parameter :: srcut = 2.0_wp*aatoau

  !> Grid size defaults
  type :: grid_defaults
    integer :: normal = gridSize(12)
    integer :: tight = gridSize(19)
    integer :: vtight = gridSize(23)
    integer :: verytight = gridSize(23)
    integer :: extreme = gridSize(32)
  end type grid_defaults
  type(grid_defaults),parameter :: lebedev = grid_defaults()

!========================================================================================!
!========================================================================================!
contains   !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!> The default setup of the surface calculator
  subroutine setup_engine(self,nat,at,xyz,model,pressure,pr,ierr,ngrid,probe,scaling,Bondi)
    !> Error source
    character(len=*),parameter :: source = 'setup_engine'

    !> Instance of the surface calculator
    class(surface_engine),intent(out) :: self

    !> Number of atoms and atom numbers
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    !> coordinates in Bohr
    real(wp),intent(in) :: xyz(3,nat)

    !> model to use
    integer, intent(in) :: model

    !> pressure
    real(wp), intent(in) :: pressure

    !> printout flag
    logical,intent(in) :: pr

    !> Error flag
    integer,intent(inout) :: ierr
    !> (optional) number of Lebedev grid points
    integer,intent(in),optional :: ngrid
    !> (optional) probe radius in Angstroem
    real(wp),intent(in),optional :: probe
    !> (optional) scaling factor for vdw radii
    real(wp),intent(in),optional :: scaling
    !> (optional) use Bondi vdW radii instead of D3
    logical,intent(in),optional :: Bondi

    !> LOCAL
    real(wp) :: probeRad
    integer :: nAng
    real(wp),allocatable :: vdwRad(:)

    !> default setting for grid size
    if (present(ngrid)) then
      nAng = ngrid
    else
      nAng = lebedev%normal
    end if

    !> default setting for probe radius
    if (present(probe)) then
      probeRad = probe*aatoau
    else
      !> 1.5 Angstroem is a typical value for the water molecule as a probe
      !> TODO: move to more accessible spot
      probeRad = 1.5_wp*aatoau
    end if
    self%probeRad_au = probeRad
    self%probeRad_aa = probeRad*autoaa

    ! TODO make not optional
    vdwRad = vanDerWaalsRadD3
    if (present(Bondi)) then
      if (Bondi) then !> Use Bondi vdW rad instead of D3?
        vdwRad = vanDerWaalsRadBondi
      end if
    end if

    if (present(scaling)) then
        vdwRad = vdwRad * scaling
    end if

    !> call init_surface_calculator.
    call self%init(at,vdwRad,probeRad,ierr,lrcut,srcut,nAng,pressure,model)

    !> call the update routine to set up neighbourlists and calculate the surface
    call self%update(at,xyz)

    if (pr) then
      call self%info(stdout)
    end if

  end subroutine setup_engine

!=========================================================================================!
!> Initialize data straucture
  subroutine init_engine(self,num,vdwRad,probeRad,ierr, &
         & rCutoff,rOffset,nAng, pressure, model)

    !> Error source
    character(len=*),parameter :: source = 'init_surface_engine'

    !> Instance of the surface calculator
    class(surface_engine),intent(inout) :: self

    !> Atomic numbers
    integer,intent(in) :: num(:)

    !> Van-der-Waals Radii
    real(wp),intent(in) :: vdwRad(:)

    !> Probe radius of the solvent
    real(wp),intent(in) :: probeRad

    !> Real-space cutoff for Born radii integration
    real(wp),intent(in) :: rCutoff

    !> Offset for surface integration cutoff
    real(wp),intent(in) :: rOffset

    !> Number of angular grid points for integration
    integer,intent(in) :: nAng

    !> pressure in a.u.
    real(wp), intent(in) :: pressure
    
    !> model for gradient calculation
    integer, intent(in) :: model

    integer :: iat,izp,jat,ij,iang
    real(wp) :: r
    integer,intent(inout) :: ierr

    self%nat = size(num)
    self%at = num
    self%ntpair = self%nat*(self%nat-1)/2
    allocate (self%ppind(2,self%ntpair))
    allocate (self%nnsas(self%nat))
    allocate (self%nnlistr(3,self%ntpair))
    allocate (self%nnlists(self%nat,self%nat))
    allocate (self%ddpair(4,self%ntpair))
    ij = 0
    do iat = 1,self%nat
      do jat = 1,iat-1
        ij = ij+1
        self%ppind(1,ij) = iat
        self%ppind(2,ij) = jat
      end do
    end do

    self%lrcut = rCutoff
    allocate (self%vdwr(self%nat))

    do iat = 1,self%nat
      izp = num(iat)
      self%vdwr(iat) = vdwRad(izp)
    end do

    allocate (self%vdwsa(self%nat))
    allocate (self%trj2(2,self%nat))
    allocate (self%wrp(self%nat))
    allocate (self%grad(3,self%nat))
    do iat = 1,self%nat
      izp = num(iat)
      self%vdwsa(iat) = vdwRad(izp)+probeRad
      self%trj2(1,iat) = (self%vdwsa(iat)-w)**2
      self%trj2(2,iat) = (self%vdwsa(iat)+w)**2
      r = self%vdwsa(iat)+w
      self%wrp(iat) = (0.25_wp/w+ &
         &            3.0_wp*ah3*(0.2_wp*r*r-0.5_wp*r*self%vdwsa(iat)+ &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
      r = self%vdwsa(iat)-w
      self%wrp(iat) = self%wrp(iat)-(0.25/w+ &
         &    3.0_wp*ah3*(0.2_wp*r*r-0.5_wp*r*self%vdwsa(iat)+ &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
    end do
    self%srcut = 2*(w+maxval(self%vdwsa))+rOffset

    self%nAng = nAng
    call bisectSearch(iang,gridSize,nAng)
    allocate (self%angGrid(3,gridSize(iang)))
    allocate (self%angWeight(gridSize(iang)))
    call getAngGrid(iang,self%angGrid,self%angWeight,ierr)

    self%model = model
    self%energy = 0.0_wp
    self%pressure = pressure

  end subroutine init_engine

!=========================================================================================!
  subroutine deallocate_engine(self)
    !> Instance of the surface calculator
    class(surface_engine),intent(inout) :: self

    if (allocated(self%ppind)) deallocate (self%ppind)
    if (allocated(self%nnsas)) deallocate (self%nnsas)
    if (allocated(self%nnlistr)) deallocate (self%nnlistr)
    if (allocated(self%nnlists)) deallocate (self%nnlists)
    if (allocated(self%ddpair)) deallocate (self%ddpair)
    if (allocated(self%vdwr)) deallocate (self%vdwr)
    if (allocated(self%vdwsa)) deallocate (self%vdwsa)
    if (allocated(self%trj2)) deallocate (self%trj2)
    if (allocated(self%wrp)) deallocate (self%wrp)
    if (allocated(self%gamsasa)) deallocate (self%gamsasa)
    if (allocated(self%grad)) deallocate (self%grad)
    if (allocated(self%angGrid)) deallocate (self%angGrid)
    if (allocated(self%angWeight)) deallocate (self%angWeight)

  end subroutine deallocate_engine

!=========================================================================================!
  subroutine update(self,num,xyz)
    !> Instance of the surface calculator
    class(surface_engine),intent(inout) :: self

    !> Atomic numbers
    integer,intent(in) :: num(:)

    !> Cartesian coordinates
    real(wp),intent(in) :: xyz(:,:)

    integer :: i

    ! initialize the neighbor list
    call update_nnlist_gbsa(self%nat,self%ntpair,self%ppind,xyz, &
       & self%lrcut,self%srcut,self%nnsas,self%nnlists,self%nnrad, &
       & self%nnlistr,self%ddpair,.false.)

      ! reset gradient und energy
       self%energy = 0.0_wp
      self%grad = 0.0_wp

    ! compute gradient and energy
    if (self%model == 0) then
      call xhcff_eg(self%nat,self%nnsas,self%nnlists,xyz,self%vdwsa, &
        & self%wrp,self%trj2,self%angWeight,self%angGrid, self%pressure, &
         & self%sasa, self%volume, self%energy, self%grad) 
    else if (self%model == 1) then
      call pv_eg(self%nat,self%nnsas,self%nnlists,xyz,self%vdwsa, &
      & self%wrp,self%trj2,self%angWeight,self%angGrid, self%pressure, &
       & self%sasa, self%volume, self%energy, self%grad) 
    end if

  end subroutine update

!=========================================================================================!
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
  subroutine update_nnlist_gbsa(nat,ntpair,ppind,xyz,lrcut,srcut, &
        & nnsas,nnlists,nnrad,nnlistr,ddpair,parallel)

    integer,intent(in) :: nat
    integer,intent(in) :: ntpair
    integer,intent(in) :: ppind(:,:)
    real(wp),intent(in) :: xyz(:,:)
    real(wp),intent(in) :: lrcut
    real(wp),intent(in) :: srcut
    integer,intent(out) :: nnsas(:)
    integer,intent(out) :: nnlists(:,:)
    integer,intent(out) :: nnrad
    integer,intent(out) :: nnlistr(:,:)
    real(wp),intent(out) :: ddpair(:,:)
    logical,intent(in) :: parallel

    if (parallel) then
      call update_nnlist_gbsa_parallel(nat,ntpair,ppind,xyz, &
         & lrcut,srcut,nnsas,nnlists,nnrad,nnlistr,ddpair)
    else
      call update_nnlist_gbsa_sequential(nat,ntpair,ppind,xyz, &
         & lrcut,srcut,nnsas,nnlists,nnrad,nnlistr,ddpair)
    end if

  end subroutine update_nnlist_gbsa

!=========================================================================================!
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
!  OMP parallel version.
  subroutine update_nnlist_gbsa_parallel(nat,ntpair,ppind,xyz,lrcut, &
        & srcut,nnsas,nnlists,nnrad,nnlistr,ddpair)
!$  use omp_lib

    integer,intent(in) :: nat
    integer,intent(in) :: ntpair
    integer,intent(in) :: ppind(:,:)
    real(wp),intent(in) :: xyz(:,:)
    real(wp),intent(in) :: lrcut
    real(wp),intent(in) :: srcut
    integer,intent(out) :: nnsas(:)
    integer,intent(out) :: nnlists(:,:)
    integer,intent(out) :: nnrad
    integer,intent(out) :: nnlistr(:,:)
    real(wp),intent(out) :: ddpair(:,:)

    integer kk,i1,i2
    real(wp) rcutn2,lrcut2,srcut2
    real(wp) x,y,z,dr2
    integer ip,ip2,thrid,nproc
    integer,allocatable :: npid(:)
    integer,allocatable :: plisttr(:,:,:)
    integer,allocatable :: nntmp(:)
    integer,allocatable :: nnls(:,:)

    nproc = 1
!$  nproc = omp_get_max_threads()

    allocate (plisttr(3,ntpair,nproc),nnls(nat,nat))
    allocate (nntmp(nat),npid(nproc))
    npid = 0

    lrcut2 = lrcut*lrcut
    srcut2 = srcut*srcut

    nnsas = 0
    nnlists = 0
!$omp parallel default(none) &
!$omp&         shared ( xyz,lrcut2,srcut2,ntpair,ppind,nat,nnlists,nnsas,ddpair ) &
!$omp&         private( i1,i2,x,y,z,dr2,ip,ip2,thrid,nntmp,nnls ) &
!$omp&         shared ( plisttr, npid )
    ip = 0
    ip2 = 0
    nntmp = 0
    nnls = 0
    thrid = 1
!$  thrid = omp_get_thread_num()+1
!$omp do
    do kk = 1,ntpair
      i1 = ppind(1,kk)
      i2 = ppind(2,kk)
      x = xyz(1,i1)-xyz(1,i2)
      y = xyz(2,i1)-xyz(2,i2)
      z = xyz(3,i1)-xyz(3,i2)
      dr2 = x**2+y**2+z**2
      ddpair(2,kk) = x
      ddpair(3,kk) = y
      ddpair(4,kk) = z
      ddpair(1,kk) = sqrt(dr2)
      if (dr2 .lt. lrcut2) then
        ip = ip+1
        plisttr(1,ip,thrid) = i1
        plisttr(2,ip,thrid) = i2
        plisttr(3,ip,thrid) = kk
        if (dr2 .lt. srcut2) then
          nntmp(i1) = nntmp(i1)+1
          nntmp(i2) = nntmp(i2)+1
          nnls(nntmp(i1),i1) = i2
          nnls(nntmp(i2),i2) = i1
        end if
      end if
    end do
!$omp end do
    npid(thrid) = ip
!$omp critical
    do i1 = 1,nat
      do i2 = 1,nntmp(i1)
        nnlists(nnsas(i1)+i2,i1) = nnls(i2,i1)
      end do
      nnsas(i1) = nnsas(i1)+nntmp(i1)
    end do
!$omp end critical
!$omp end parallel

    nnrad = 0
    do thrid = 1,nproc
      do kk = nnrad+1,nnrad+npid(thrid)
        nnlistr(1:3,kk) = plisttr(1:3,kk-nnrad,thrid)
      end do
      nnrad = nnrad+npid(thrid)
    end do

    deallocate (nntmp,nnls)

  end subroutine update_nnlist_gbsa_parallel

!=========================================================================================!
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut
!  Sequential version.
  pure subroutine update_nnlist_gbsa_sequential(nat,ntpair,ppind,xyz,lrcut, &
        & srcut,nnsas,nnlists,nnrad,nnlistr,ddpair)

    integer,intent(in) :: nat
    integer,intent(in) :: ntpair
    integer,intent(in) :: ppind(:,:)
    real(wp),intent(in) :: xyz(:,:)
    real(wp),intent(in) :: lrcut
    real(wp),intent(in) :: srcut
    integer,intent(out) :: nnsas(:)
    integer,intent(out) :: nnlists(:,:)
    integer,intent(out) :: nnrad
    integer,intent(out) :: nnlistr(:,:)
    real(wp),intent(out) :: ddpair(:,:)

    integer kk,i1,i2
    real(wp) rcutn2,lrcut2,srcut2
    real(wp) x,y,z,dr2
    integer ip,ip2

    lrcut2 = lrcut*lrcut
    srcut2 = srcut*srcut

    nnsas = 0
    nnlists = 0
    ip = 0
    ip2 = 0
    nnlistr = 0
    do kk = 1,ntpair
      i1 = ppind(1,kk)
      i2 = ppind(2,kk)
      x = xyz(1,i1)-xyz(1,i2)
      y = xyz(2,i1)-xyz(2,i2)
      z = xyz(3,i1)-xyz(3,i2)
      dr2 = x**2+y**2+z**2
      ddpair(2,kk) = x
      ddpair(3,kk) = y
      ddpair(4,kk) = z
      ddpair(1,kk) = sqrt(dr2)
      if (dr2 .lt. lrcut2) then
        ip = ip+1
        nnlistr(1,ip) = i1
        nnlistr(2,ip) = i2
        nnlistr(3,ip) = kk
        if (dr2 .lt. srcut2) then
          nnsas(i1) = nnsas(i1)+1
          nnsas(i2) = nnsas(i2)+1
          nnlists(nnsas(i1),i1) = i2
          nnlists(nnsas(i2),i2) = i1
        end if
      end if
    end do
    nnrad = ip

  end subroutine update_nnlist_gbsa_sequential

!=========================================================================================!
  subroutine info(self,unit)
    !> Data structure
    class(surface_engine),intent(in) :: self

    !> Unit for IO
    integer,intent(in) :: unit

   ! if (allocated(self%gamsasa)) then
   !   write (unit,'(2x, a, t40, es14.4, 1x, a, t60, es14.4, 1x, a)') &
   !     "Surface tension",surfaceTension,"Eh",surfaceTension*automNm,"dyn/cm"
   ! end if

    write (unit,'(2x, a, t40, i14, 1x, a)') &
      "Grid points",self%nAng,"/ per atom"

    write (unit,'(2x, a, t40, f14.4, 1x, a)') &
      "Probe radius",self%probeRad_au,"/ Bohr"

    write (unit,'(2x, a )') repeat('-',54)

 !   Todo: use area variable
    write (unit,'(2x, a, t40, f14.4, 1x, a)') &
      "Total molecular surface",self%sasa,"/ Bohr²"
    write (unit,'(2x, a, t40, f14.4, 1x, a)') &
      "Total molecular surface",self%sasa*autoaa**2,"/ Ang²"

    write (unit,*)
  end subroutine info

!=========================================================================================!
!=========================================================================================!
end module libpv_surface_engine
