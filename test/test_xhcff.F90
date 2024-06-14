module test_xhcff
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use xhcfflib_interface
  use xhcff_surface_module
  use xhcff_type_timer
  implicit none
  private

  public :: collect_xhcff

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 1e-6

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for XHCFF calculations
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_xhcff(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("XHCFF singlepoint          ",test_xhcff_sp), &
    new_unittest("XHCFF (Bondi radii) SP     ",test_xhcff_bondi_sp), &
    new_unittest("XHCFF OpenMP parallel SP   ",test_xhcff_openmp), &
    new_unittest("XHCFF(+PV) singlepoint     ",test_xhcff_pv_sp), &
    new_unittest("XHCFF(+PV) num. gradient   ",test_xhcff_pv_numgrad) &
    ]
!&>
  end subroutine collect_xhcff

!========================================================================================!

  subroutine test_xhcff_sp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    real(wp) :: p,probe
    type(xhcfflib_calculator) :: xhcff
!&<
    real(wp),parameter :: e_ref = 0.0_wp  !> no energy for XHCFF standard implementation
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    &  -0.001050219817793_wp,   0.000044408915900_wp,  -0.000000111886661_wp, &
    &  -0.000395314094809_wp,   0.000023296469152_wp,  -0.000003586648768_wp, &
    &  -0.000464946228945_wp,   0.002016143817697_wp,  -0.000000398820235_wp, &
    &   0.004579861978900_wp,   0.003229979134347_wp,  -0.000001207175299_wp, &
    &   0.001186382382069_wp,   0.000201440040887_wp,  -0.000002036876930_wp, &
    &  -0.001114973497024_wp,  -0.000277940089010_wp,  -0.000004536452373_wp, &
    &  -0.000691391833646_wp,  -0.000589888711945_wp,  -0.000003944601526_wp, &
    &  -0.005080163909164_wp,  -0.001648839344471_wp,  -0.000010075220713_wp, &
    &  -0.000024110958398_wp,   0.000087929161037_wp,   0.000000583890724_wp, &
    &   0.000397717516352_wp,  -0.000225005149492_wp,  -0.000002475952831_wp, &
    &   0.003994304635340_wp,  -0.003537552592492_wp,  -0.000000105536085_wp, &
    &   0.000122037564703_wp,  -0.000001583053777_wp,  -0.000006032398460_wp, &
    &   0.000945738206116_wp,   0.000239975834993_wp,  -0.000001384068474_wp, &
    &  -0.000531458005767_wp,  -0.001008059492273_wp,   0.000003864347621_wp, &
    &  -0.001423126810413_wp,   0.004222474230492_wp,   0.000002825628412_wp, &
    &  -0.001493343440606_wp,  -0.002056960897210_wp,   0.003655854266956_wp, &
    &  -0.001493326508727_wp,  -0.002052108033756_wp,  -0.003658677775602_wp, &
    &  -0.001475218162121_wp,   0.004194096903855_wp,   0.000002248102452_wp, &
    &   0.003189000054599_wp,  -0.003063135710907_wp,  -0.000003097494498_wp, &
    &   0.000573937255807_wp,   0.002505385805384_wp,   0.003634723790020_wp, &
    &   0.000571424492870_wp,   0.002509701193112_wp,  -0.003632244425056_wp, &
    &  -0.004210743378849_wp,  -0.001347960727226_wp,   0.000014758407531_wp, &
    &   0.001948244596208_wp,  -0.001712525271071_wp,   0.003616726552258_wp, &
    &   0.001939687963298_wp,  -0.001753272433224_wp,  -0.003601669652463_wp &
    & ], shape(g_ref))
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    p = 10.0_wp     !> pressure in GPa
    probe = 0.0_wp  !> Probe radius
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

    !> calculation
    call xhcff%init(nat,at,xyz,p,'XHCFF',proberad=probe,verbose=.true.,printlevel=2)
    call xhcff%singlepoint(nat,at,xyz,energy,grad,iostat=io)
    !write (*,'(F25.15)') energy
    !write (*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-7_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient does not match reference")
     call xhcff%info()
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_xhcff_sp

!=======================================================================================!

  subroutine test_xhcff_bondi_sp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    real(wp) :: p,probe
    type(xhcfflib_calculator) :: xhcff
!&<
    real(wp),parameter :: e_ref = 0.0_wp  !> no energy for XHCFF standard implementation
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    &  -0.004413376075752_wp,   0.000954403196482_wp,  -0.000000161399853_wp, &
    &  -0.000258787756931_wp,   0.000031462790573_wp,  -0.000004387443464_wp, &
    &  -0.000816843814755_wp,   0.004555283033469_wp,   0.000000560666814_wp, &
    &   0.005781419585546_wp,   0.004561759669946_wp,  -0.000001562830498_wp, &
    &   0.000942099034853_wp,   0.000205292230599_wp,  -0.000002891292393_wp, &
    &  -0.000812150631670_wp,  -0.000118070779719_wp,  -0.000005944253768_wp, &
    &  -0.000464117954306_wp,  -0.000471953680629_wp,  -0.000005367852321_wp, &
    &  -0.007874676690601_wp,  -0.002655555289088_wp,  -0.000019549878554_wp, &
    &  -0.000012951965609_wp,   0.000039545605676_wp,   0.000001277457777_wp, &
    &   0.000219046865019_wp,  -0.000076955857305_wp,  -0.000005261393871_wp, &
    &   0.006325604624626_wp,  -0.005456670291075_wp,  -0.000002771304701_wp, &
    &   0.000084492940621_wp,   0.000002994123469_wp,  -0.000007186240015_wp, &
    &   0.004225682470689_wp,   0.001333386588531_wp,  -0.000003504364648_wp, &
    &  -0.001463840915938_wp,  -0.004367415920492_wp,   0.000023776384410_wp, &
    &  -0.001406105001708_wp,   0.004068631101217_wp,   0.000002771680378_wp, &
    &  -0.001443870801096_wp,  -0.001984302372547_wp,   0.003539190491154_wp, &
    &  -0.001443792167358_wp,  -0.001979736313148_wp,  -0.003541777484641_wp, &
    &  -0.001421271852337_wp,   0.004045355769066_wp,   0.000002200487454_wp, &
    &   0.003193141894393_wp,  -0.002586873593383_wp,  -0.000003313948886_wp, &
    &   0.000565464482660_wp,   0.002428162148809_wp,   0.003521095348850_wp, &
    &   0.000563062811016_wp,   0.002432305506726_wp,  -0.003518555845146_wp, &
    &  -0.003789954827387_wp,  -0.001568821045609_wp,   0.000016365907209_wp, &
    &   0.001865168722117_wp,  -0.001676305102592_wp,   0.003503075221728_wp, &
    &   0.001856557023907_wp,  -0.001715921518977_wp,  -0.003488078113014_wp &
    & ], shape(g_ref))
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    p = 10.0_wp     !> pressure in GPa
    probe = 0.0_wp  !> Probe radius
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

    !> calculation
    call xhcff%init(nat,at,xyz,p,'XHCFF',proberad=probe,gridpts=5294,vdwSet=1, &
    &    verbose=.false.,printlevel=2)
    call xhcff%singlepoint(nat,at,xyz,energy,grad,iostat=io)
    !write (*,'(F25.15)') energy
    !write (*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-7_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient does not match reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_xhcff_bondi_sp

!========================================================================================!

  subroutine test_xhcff_openmp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,i,j,ntimes
    real(wp) :: p,probe
    type(xhcfflib_calculator) :: xhcff
    type(xhcff_timer) :: timer
    character(len=40) :: atmp
    logical :: speedup
!&<
    real(wp),parameter :: e_ref = 0.0_wp  !> no energy for XHCFF standard implementation
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    &  -0.004413376075752_wp,   0.000954403196482_wp,  -0.000000161399853_wp, &
    &  -0.000258787756931_wp,   0.000031462790573_wp,  -0.000004387443464_wp, &
    &  -0.000816843814755_wp,   0.004555283033469_wp,   0.000000560666814_wp, &
    &   0.005781419585546_wp,   0.004561759669946_wp,  -0.000001562830498_wp, &
    &   0.000942099034853_wp,   0.000205292230599_wp,  -0.000002891292393_wp, &
    &  -0.000812150631670_wp,  -0.000118070779719_wp,  -0.000005944253768_wp, &
    &  -0.000464117954306_wp,  -0.000471953680629_wp,  -0.000005367852321_wp, &
    &  -0.007874676690601_wp,  -0.002655555289088_wp,  -0.000019549878554_wp, &
    &  -0.000012951965609_wp,   0.000039545605676_wp,   0.000001277457777_wp, &
    &   0.000219046865019_wp,  -0.000076955857305_wp,  -0.000005261393871_wp, &
    &   0.006325604624626_wp,  -0.005456670291075_wp,  -0.000002771304701_wp, &
    &   0.000084492940621_wp,   0.000002994123469_wp,  -0.000007186240015_wp, &
    &   0.004225682470689_wp,   0.001333386588531_wp,  -0.000003504364648_wp, &
    &  -0.001463840915938_wp,  -0.004367415920492_wp,   0.000023776384410_wp, &
    &  -0.001406105001708_wp,   0.004068631101217_wp,   0.000002771680378_wp, &
    &  -0.001443870801096_wp,  -0.001984302372547_wp,   0.003539190491154_wp, &
    &  -0.001443792167358_wp,  -0.001979736313148_wp,  -0.003541777484641_wp, &
    &  -0.001421271852337_wp,   0.004045355769066_wp,   0.000002200487454_wp, &
    &   0.003193141894393_wp,  -0.002586873593383_wp,  -0.000003313948886_wp, &
    &   0.000565464482660_wp,   0.002428162148809_wp,   0.003521095348850_wp, &
    &   0.000563062811016_wp,   0.002432305506726_wp,  -0.003518555845146_wp, &
    &  -0.003789954827387_wp,  -0.001568821045609_wp,   0.000016365907209_wp, &
    &   0.001865168722117_wp,  -0.001676305102592_wp,   0.003503075221728_wp, &
    &   0.001856557023907_wp,  -0.001715921518977_wp,  -0.003488078113014_wp &
    & ], shape(g_ref))
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    p = 10.0_wp     !> pressure in GPa
    probe = 0.0_wp  !> Probe radius
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

#ifdef WITH_OpenMP

    ntimes = 2
    call timer%new(ntimes,.true.)
    do i = 1,ntimes
      call OMP_Set_Num_Threads(i)
#ifdef WITH_MKL
      call MKL_Set_Num_Threads(i)
#endif
      call ompprint_intern(atmp)

      call timer%measure(i,atmp)

      !> calculation
      call xhcff%reset
      grad = 0.0_wp
      do j = 1,20
        call xhcff%init(nat,at,xyz,p,'XHCFF',proberad=probe,gridpts=5294,vdwSet=1, &
        &    verbose=.false.,printlevel=2)
        call xhcff%singlepoint(nat,at,xyz,energy,grad,iostat=io)
        !write (*,'(F25.15)') energy
        !write (*,'(3(F20.15,"_wp,")," &")') grad
        call check(error,io,0)
        if (allocated(error)) return

        call check(error,energy,e_ref,thr=1e-7_wp)
        if (allocated(error)) return

        if (any(abs(grad-g_ref) > thr2)) then
          call test_failed(error,"Gradient does not match reference")
          print'(3es21.14)',grad
          print'("---")'
          print'(3es21.14)',g_ref
          print'("---")'
          print'(3es21.14)',grad-g_ref
        end if
      end do

      call timer%measure(i)
      !write (*,*)
      !call timer%write_timing(stdout,i)
    end do

    call OMP_Set_Num_Threads(1)
#ifdef WITH_MKL
    call MKL_Set_Num_Threads(1)
#endif

    !> check for any noticable speedup
    speedup = 1.1_wp < (timer%get(1)/timer%get(2))
#else /* WITH_OpenMP */
    speedup = .true.
#endif /* WITH_OpenMP */
    call check(error,speedup,.true.)
    if (allocated(error)) return

    deallocate (grad)
  end subroutine test_xhcff_openmp

  subroutine ompprint_intern(str)
    use omp_lib
    implicit none
    integer :: nproc,TID
    character(len=*) :: str
!$OMP PARALLEL PRIVATE(TID)
    TID = OMP_GET_THREAD_NUM()
    IF (TID .EQ. 0) THEN
      nproc = OMP_GET_NUM_THREADS()
#ifdef WITH_MKL
      write (str,'(a,i3)') 'OMP/MKL threads = ',nproc
#else
      write (str,'(a,i3)') 'OMP threads = ',nproc
#endif
    END IF
!$OMP END PARALLEL
  end subroutine ompprint_intern

!========================================================================================!

  subroutine test_xhcff_pv_sp(error)
    use coffeine
    use pv_interface
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    real(wp) :: p,probe
    type(xhcfflib_calculator) :: pv

!&<
    real(wp),parameter :: e_ref = 1053.433399736551337_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    & -12.524435266288961_wp,   2.825689169961084_wp,   0.007355745321851_wp, &
    &  -0.548359748562616_wp,  -0.461096971187015_wp,   0.009651060024959_wp, &
    &  -2.946963997140400_wp,  13.461943476494234_wp,  -0.020885848627251_wp, &
    &  17.477529661105546_wp,  13.907195035829664_wp,   0.002926844227702_wp, &
    &   2.534481924561555_wp,   0.427557269647703_wp,  -0.015761088339607_wp, &
    &  -3.012198583232461_wp,  -0.514333513985023_wp,  -0.021965043017575_wp, &
    &  -1.069657399318409_wp,  -2.442757763760082_wp,  -0.021256584022544_wp, &
    & -23.704445258929475_wp,  -6.624818119146281_wp,  -0.051421817833663_wp, &
    &  -0.121023125688725_wp,   0.502170376262657_wp,   0.030822831277590_wp, &
    &   0.258606657848168_wp,   0.165682877793313_wp,  -0.028194769886269_wp, &
    &  18.807113353287065_wp, -17.911396655425676_wp,  -0.007297948460879_wp, &
    &  -0.354228923875067_wp,   0.235846316471151_wp,  -0.000115613244500_wp, &
    &  12.436798867600693_wp,   5.019142891186579_wp,   0.003201449604061_wp, &
    &  -5.896292596214175_wp, -13.080985912928634_wp,   0.088690786797322_wp, &
    &  -3.938088193296592_wp,  13.078900430379376_wp,   0.019305469534641_wp, &
    &  -4.396663972059278_wp,  -6.731798238601612_wp,  10.987364041575130_wp, &
    &  -4.426315050406996_wp,  -6.744838716823117_wp, -11.029560194508701_wp, &
    &  -4.916084877899209_wp,  12.754280865805390_wp,   0.017661353744272_wp, &
    &  11.090270805982300_wp,  -9.889077434457262_wp,  -0.002372703129355_wp, &
    &   1.412520560538434_wp,   7.731695501777458_wp,  11.456904354301424_wp, &
    &   1.420835074492471_wp,   7.780791085679365_wp, -11.480929088914742_wp, &
    & -12.134804883456233_wp,  -3.942805670782521_wp,   0.061630641745890_wp, &
    &   7.295746804481155_wp,  -4.734743086698068_wp,  12.545592707503365_wp, &
    &   7.308897091683217_wp,  -4.909208690446083_wp, -12.551943787039395_wp &
    & ], shape(g_ref))
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    p =  1.0_wp/3.4e-5_wp    !> set pressure to 1 a.u., to directly compare volume gradient
    probe = 0.0_wp  !> Probe radius
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)

    !> calculation
    call pv%init(nat,at,xyz,p,'PV',gridpts=5294, &
    &    proberad=probe,vdwSet=1,verbose=.false.,printlevel=2)
    call pv%singlepoint(nat,at,xyz,energy,grad,iostat=io)
    !call pv%info()
    !write (*,'(F25.15)') energy
    !write (*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=5e-4_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient does not match reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_xhcff_pv_sp

!========================================================================================!

  subroutine test_xhcff_pv_numgrad(error)
    use coffeine
    use pv_interface
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,i,j
    real(wp) :: p,probe,step,bw,bw2,fw,fw2
    type(xhcfflib_calculator) :: pv
    real(wp),allocatable :: gradient(:,:),g_ref(:,:),stencil(:,:)
!&<
    real(wp),parameter :: e_ref = 1053.415994538545_wp
!&>

    !> setup
    nat = testnat
    allocate (at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    p =  1.0_wp/3.4e-5_wp    !> set pressure to 1 a.u., to directly compare volume gradient
    probe = 0.0_wp  !> Probe radius
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)
    allocate (gradient(3,nat),g_ref(3,nat),stencil(3,nat),source=0.0_wp)

    !> calculation
    call pv%init(nat,at,xyz,p,'PV',gridpts=974, &
    &    proberad=probe,vdwSet=1,verbose=.false.,printlevel=2)
    call pv%singlepoint(nat,at,xyz,energy,grad,iostat=io)
    !call pv%info()
    !write (*,'(F25.15)') energy
    !write (*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    stencil = xyz
    step = 0.001_wp
    do i = 1,nat
      do j = 1,3
        !write (*,*) 'Numerical gradient dimension ', (i-1)*3+j
        stencil(j,i) = stencil(j,i)-2.0_wp*step
        call pv%singlepoint(nat,at,stencil,bw2,gradient)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)-1.0_wp*step
        call pv%singlepoint(nat,at,stencil,bw,gradient)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)+1.0_wp*step
        call pv%singlepoint(nat,at,stencil,fw,gradient)
        stencil(j,i) = xyz(j,i)
        stencil(j,i) = stencil(j,i)+2.0_wp*step
        call pv%singlepoint(nat,at,stencil,fw2,gradient)
        stencil(j,i) = xyz(j,i)
        g_ref(j,i) = (bw2/12.0_wp-8.0_wp*bw/12.0_wp+8.0_wp*fw/12.0_wp-fw2/12.0_wp)/step
      end do
    end do

    call check(error,energy,e_ref,thr=5e-4_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > 5e-4_wp)) then
      call test_failed(error,"Gradient does not match reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_xhcff_pv_numgrad

!========================================================================================!
!========================================================================================!
end module test_xhcff
