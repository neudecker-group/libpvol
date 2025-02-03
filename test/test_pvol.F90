module test_pvol
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use libpvol_interface
  use pvol_type_timer
  implicit none
  private

  public :: collect_pvol

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 1e-6

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for PV calculations
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_pvol(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("XHCFF(+PV) singlepoint          ",test_xhcff_sp), &
    new_unittest("XHCFF(+PV) (Bondi radii) SP     ",test_xhcff_bondi_sp), &
    new_unittest("XHCFF+(PV) OpenMP parallel SP   ",test_xhcff_openmp), &
    new_unittest("PV singlepoint     ",test_pv_sp), &
    new_unittest("PV num. gradient   ",test_pv_numgrad), &
    new_unittest("probe and scaling variables  ",test_probe_scaling) &
    ]
!&>
  end subroutine collect_pvol

!========================================================================================!

  subroutine test_xhcff_sp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    real(wp) :: p,probe
    type(libpvol_calculator) :: xhcff
!&<
    real(wp),parameter :: e_ref = 0.28384347661327025_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    & -0.001041660129153_wp,   0.000044012240474_wp,  -0.000000104495569_wp, &
    & -0.000391515818108_wp,   0.000023035254661_wp,  -0.000003543125664_wp, &
    & -0.000461319864782_wp,   0.001999127208894_wp,  -0.000000389010553_wp, &
    &  0.004532115201258_wp,   0.003196495249987_wp,  -0.000001188304606_wp, &
    &  0.001176089546712_wp,   0.000199719735782_wp,  -0.000002013259920_wp, &
    & -0.001105867993714_wp,  -0.000275619655264_wp,  -0.000004491766213_wp, &
    & -0.000685856738822_wp,  -0.000584938834796_wp,  -0.000003904904131_wp, &
    & -0.005021760639710_wp,  -0.001629804598888_wp,  -0.000009952485433_wp, &
    & -0.000024158512589_wp,   0.000086998322885_wp,   0.000000584205170_wp, &
    &  0.000394072409102_wp,  -0.000223130909181_wp,  -0.000002448634810_wp, &
    &  0.003947848484389_wp,  -0.003496692946031_wp,  -0.000000098004834_wp, &
    &  0.000120475846298_wp,  -0.000001586505257_wp,  -0.000005963536293_wp, &
    &  0.000937473782473_wp,   0.000237930708164_wp,  -0.000001365954045_wp, &
    & -0.000527271007769_wp,  -0.000999584815465_wp,   0.000003838222672_wp, &
    & -0.001402255480427_wp,   0.004159619622558_wp,   0.000002789742984_wp, &
    & -0.001471427117831_wp,  -0.002026362367241_wp,   0.003601452505907_wp, &
    & -0.001471410437940_wp,  -0.002021581726098_wp,  -0.003604221660462_wp, &
    & -0.001453571587200_wp,   0.004131664618975_wp,   0.000002220812003_wp, &
    &  0.003141231814437_wp,  -0.003017562873264_wp,  -0.000003045229573_wp, &
    &  0.000565087456459_wp,   0.002468085614068_wp,   0.003580636501425_wp, &
    &  0.000562612089493_wp,   0.002472336778421_wp,  -0.003578181701917_wp, &
    & -0.004148385591833_wp,  -0.001327913829580_wp,   0.000014544933573_wp, &
    &  0.001918941789417_wp,  -0.001687052773844_wp,   0.003562907105962_wp, &
    &  0.001910512499839_wp,  -0.001727193519961_wp,  -0.003548061955673_wp  &
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
    call xhcff%init(nat,at,xyz,p,'XHCFF',proberad=probe,verbose=.false.,printlevel=2)
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
    type(libpvol_calculator) :: xhcff
!&<
    real(wp),parameter :: e_ref = 0.354862618399131_wp  
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    & -0.004386373900270_wp,   0.000948715704419_wp,  -0.000000150913485_wp, &
    & -0.000257181940933_wp,   0.000031454048740_wp,  -0.000004345297355_wp, &
    & -0.000812103530847_wp,   0.004527306739113_wp,   0.000000566683718_wp, &
    &  0.005738107520046_wp,   0.004528061465676_wp,  -0.000001541689568_wp, &
    &  0.000935951777051_wp,   0.000204241601461_wp,  -0.000002863908479_wp, &
    & -0.000807439397726_wp,  -0.000117119849646_wp,  -0.000005897972569_wp, &
    & -0.000461560980401_wp,  -0.000468812278680_wp,  -0.000005325138943_wp, &
    & -0.007814119794548_wp,  -0.002634799729654_wp,  -0.000019389219228_wp, &
    & -0.000013174302573_wp,   0.000039476756149_wp,   0.000001277476960_wp, &
    &  0.000217375170854_wp,  -0.000076259421250_wp,  -0.000005219339453_wp, &
    &  0.006276384393177_wp,  -0.005414260346640_wp,  -0.000002740354330_wp, &
    &  0.000083545956555_wp,   0.000003197088311_wp,  -0.000007123280712_wp, &
    &  0.004199210418579_wp,   0.001325353257039_wp,  -0.000003473185927_wp, &
    & -0.001455095835668_wp,  -0.004340162303853_wp,   0.000023638700173_wp, &
    & -0.001385829335831_wp,   0.004009248176674_wp,   0.000002740755404_wp, &
    & -0.001423041564869_wp,  -0.001954961516530_wp,   0.003487322808366_wp, &
    & -0.001422964083743_wp,  -0.001950462386340_wp,  -0.003489852477230_wp, &
    & -0.001400773871060_wp,   0.003986314013886_wp,   0.000002177935014_wp, &
    &  0.003146001889155_wp,  -0.002548700260187_wp,  -0.000003255670826_wp, &
    &  0.000556840921568_wp,   0.002392825186109_wp,   0.003469492904304_wp, &
    &  0.000554474453579_wp,   0.002396907810767_wp,  -0.003466971220072_wp, &
    & -0.003734736737355_wp,  -0.001545570306762_wp,   0.000016135718323_wp, &
    &  0.001837494121684_wp,  -0.001651478864538_wp,   0.003451736915854_wp, &
    &  0.001829008653578_wp,  -0.001690514584266_wp,  -0.003436940229942_wp  &
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
    type(libpvol_calculator) :: xhcff
    type(pvol_timer) :: timer
    character(len=40) :: atmp
    logical :: speedup
!&<
    real(wp),parameter :: e_ref = 0.354862618399131_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    & -0.004386373900270_wp,   0.000948715704419_wp,  -0.000000150913485_wp, &
    & -0.000257181940933_wp,   0.000031454048740_wp,  -0.000004345297355_wp, &
    & -0.000812103530847_wp,   0.004527306739113_wp,   0.000000566683718_wp, &
    &  0.005738107520046_wp,   0.004528061465676_wp,  -0.000001541689568_wp, &
    &  0.000935951777051_wp,   0.000204241601461_wp,  -0.000002863908479_wp, &
    & -0.000807439397726_wp,  -0.000117119849646_wp,  -0.000005897972569_wp, &
    & -0.000461560980401_wp,  -0.000468812278680_wp,  -0.000005325138943_wp, &
    & -0.007814119794548_wp,  -0.002634799729654_wp,  -0.000019389219228_wp, &
    & -0.000013174302573_wp,   0.000039476756149_wp,   0.000001277476960_wp, &
    &  0.000217375170854_wp,  -0.000076259421250_wp,  -0.000005219339453_wp, &
    &  0.006276384393177_wp,  -0.005414260346640_wp,  -0.000002740354330_wp, &
    &  0.000083545956555_wp,   0.000003197088311_wp,  -0.000007123280712_wp, &
    &  0.004199210418579_wp,   0.001325353257039_wp,  -0.000003473185927_wp, &
    & -0.001455095835668_wp,  -0.004340162303853_wp,   0.000023638700173_wp, &
    & -0.001385829335831_wp,   0.004009248176674_wp,   0.000002740755404_wp, &
    & -0.001423041564869_wp,  -0.001954961516530_wp,   0.003487322808366_wp, &
    & -0.001422964083743_wp,  -0.001950462386340_wp,  -0.003489852477230_wp, &
    & -0.001400773871060_wp,   0.003986314013886_wp,   0.000002177935014_wp, &
    &  0.003146001889155_wp,  -0.002548700260187_wp,  -0.000003255670826_wp, &
    &  0.000556840921568_wp,   0.002392825186109_wp,   0.003469492904304_wp, &
    &  0.000554474453579_wp,   0.002396907810767_wp,  -0.003466971220072_wp, &
    & -0.003734736737355_wp,  -0.001545570306762_wp,   0.000016135718323_wp, &
    &  0.001837494121684_wp,  -0.001651478864538_wp,   0.003451736915854_wp, &
    &  0.001829008653578_wp,  -0.001690514584266_wp,  -0.003436940229942_wp  &
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

  subroutine test_pv_sp(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io
    real(wp) :: p,probe
    type(libpvol_calculator) :: pv

!&<
    real(wp),parameter :: e_ref = 1043.713583526855246_wp
    real(wp),parameter :: g_ref(3,testnat) = reshape([&
    & -12.465606743005154_wp,   2.817608030716889_wp,   0.000816706280871_wp, &
    & -0.543891620880669_wp,  -0.460105502956387_wp,   0.010676250684199_wp, &
    & -2.862862832889010_wp,  13.532214959882120_wp,  -0.023996430281781_wp, &
    & 17.348119994783627_wp,  13.783744709848785_wp,   0.003775494983762_wp, &
    &  2.526856228441928_wp,   0.437327004284599_wp,  -0.016434820148027_wp, &
    & -2.996950561223764_wp,  -0.505875992397041_wp,  -0.022232147210393_wp, &
    & -1.069800936959756_wp,  -2.439561709892783_wp,  -0.022442690259254_wp, &
    & -23.536070656472845_wp,  -6.512801480503610_wp,  -0.050756209053277_wp, &
    & -0.117749539506414_wp,   0.498585993231667_wp,   0.031805180579671_wp, &
    &  0.258022212019399_wp,   0.166871995614419_wp,  -0.029451851210697_wp, &
    & 18.685599417429891_wp, -17.778415205085274_wp,  -0.006677761231714_wp, &
    & -0.353673512551911_wp,   0.226719724993002_wp,   0.000982649441637_wp, &
    & 13.077147500184282_wp,   4.968686766552604_wp,  -0.003225050209747_wp, &
    & -5.628174432985797_wp, -13.516481382850161_wp,   0.083635570326213_wp, &
    & -3.884989977834100_wp,  12.854947572412016_wp,   0.021562034244682_wp, &
    & -4.348048994076787_wp,  -6.620087726323123_wp,  10.793023571145223_wp, &
    & -4.377934507479155_wp,  -6.633891403452867_wp, -10.831184242245213_wp, &
    & -4.943450119588642_wp,  12.493344118081943_wp,   0.019946953519762_wp, &
    & 10.715659013513942_wp,  -9.733588497589315_wp,   0.000206114925098_wp, &
    &  1.167752708422402_wp,   7.631634206692411_wp,  11.265032447910441_wp, &
    &  1.176259790015803_wp,   7.681048506408681_wp, -11.285287690168893_wp, &
    & -11.941574514980100_wp,  -3.806874616913637_wp,   0.063427689661481_wp, &
    &  7.080788982857712_wp,  -4.506885759430385_wp,  12.378845423064384_wp, &
    &  7.094913858903276_wp,  -4.679708249541258_wp, -12.382863930870178_wp & 
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
  end subroutine test_pv_sp

!========================================================================================!

  subroutine test_pv_numgrad(error)
    use coffeine
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: energy
    real(wp),allocatable :: xyz(:,:),grad(:,:)
    integer,allocatable :: at(:)
    integer :: nat,io,i,j
    real(wp) :: p,probe,step,bw,bw2,fw,fw2
    type(libpvol_calculator) :: pv
    real(wp),allocatable :: gradient(:,:),g_ref(:,:),stencil(:,:)
!&<
    real(wp),parameter :: e_ref = 1043.6963300529367_wp
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
   ! call pv%info()
   ! write (*,'(F25.15)') energy
   ! write (*,'(3(F20.15,"_wp,")," &")') grad
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
  end subroutine test_pv_numgrad

  !> also tests if code is close to analytical volume of a sphere!
  subroutine test_probe_scaling(error)  
    implicit none
    type(error_type),allocatable,intent(out) :: error
    !> testparameters
    integer :: nat = 1
    integer :: at(1) = [2]

    real(wp), allocatable :: xyz(:,:)
    real(wp) :: p = 1.0_wp/3.3989309735473356e-5_wp !> set pressure to 1 a.u., to directly compare volume
    real(wp) :: probe = 1.5 !> testproberad in Angstrom
    real(wp) :: scaling = 1.5 !> testscale 
    real(wp) :: energy
    integer :: io
    type(libpvol_calculator) :: pv
    real(wp),allocatable :: grad(:,:)
    real(wp) :: refvol_probe = 689.4123885
    real(wp) :: refvol_scale = 261.78392243

    !> setup
    energy = 0.0_wp
    allocate (grad(3,nat),source=0.0_wp)
    allocate (xyz(3,nat),source=0.0_wp)

    !> calculation with proberadius
    call pv%init(nat,at,xyz,p,'PV',gridpts=5294, &
    &    proberad=probe,vdwSet=1,verbose=.false.,printlevel=2)
    call pv%singlepoint(nat,at,xyz,energy,grad,iostat=io)

    call check(error,io,0)
    if (allocated(error)) return

    !> at chosen pressure energy is equal to volume
    call check(error,energy,refvol_probe,thr=5e-4_wp)
    if (allocated(error)) return

    call pv%reset()
    !> calculation with scaling
    call pv%init(nat,at,xyz,p,'PV',gridpts=5294, &
    &    proberad=0.0_wp, scaling=scaling, vdwSet=1,verbose=.false.,printlevel=2)
    call pv%singlepoint(nat,at,xyz,energy,grad,iostat=io)

    call check(error,io,0)
    if (allocated(error)) return

    !> at chosen pressure energy is equal to volume
    call check(error,energy,refvol_scale,thr=5e-4_wp)
    if (allocated(error)) return
  end subroutine
!========================================================================================!
end module test_pvol
