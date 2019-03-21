program VORONOI

  !====================
  ! To compile and run
  ! make
  ! ./freel
  !===================

  use mesh_data
  use alloc
  use bibli_init
  use our_module

  implicit none
  !
  ! Declarations
  ! ------------
  ! Arrays
  real*8,  dimension(:,:),   allocatable :: XYp,XYp0,generator_length,XY_c, Weight,grad
  real*8,  dimension(:),     allocatable :: Vp, V_c, V_cn, V_n, W_c,WW_c,W_p
  integer, dimension(:,:),   allocatable :: Cp
  integer, dimension(:),     allocatable :: Npa
  ! Variables
  integer  :: ncells, npoints, nx, ny, Nppc, npart, npart0
  integer  :: nstep,channel,exxit,problem,method
  ! Characters
  character(10)       :: filename,filename2,filename3,filename33
  character(13)       :: stringg,filename11,filename22
  character(30)       :: cnum,filefct,cnum2
  !-----------
  ! Reals
  real*8  :: critical_length,critical_angle
  real*8  :: S1,S2,GG,dl,xx,yy,energy,energy0, pot

  ! Integer
  integer :: i,j,k,kk,jj,ii,in,inn,ne, p
  integer :: maxcycle, max_points,max_number_generators
  integer, external :: iDirGen
  !
  print*,'************************************************'
  print*,'*                                              *'
  print*,'*   V O R O N O I   M E S H E R    2 D         *'
  print*,'*                                              *'
  print*,'************************************************'
  !
  ! Debug=0, quiet. Debug=1, print to screen
  debug = 0

  !----------------------------------------------------
  ! Define small real numbers to compare two reals
  eps = epsilon(eps) * 100d0
  eps = 1.e-10
  !----------------------------------------------------


  !=========================================================
  ! Set the problem number to run (see bibli_init.f90 file)
  print*,' Which initialisation of generators? (0:random,1:uniform,-1:from file)';
  read*, problem
  !if( abs(problem)>1 ) problem = 0
  !=========================================================

  ! Computational domain
  !print*,' Domain to test: '
  j = 0
  select case( j )
  case(2)
    print*,' Enter xmin, xmax, ymin, ymax'
    read*,xmin,xmax,ymin,ymax
  case(1)
    xmin = -1.1_rr ; ymin = -1.2_rr
    xmax =  6.0_rr ; ymax =  4.5_rr
  case(0)
    xmin = 0.0_rr ; ymin = 0.0_rr
    xmax = 1.0_rr ; ymax = 1.0_rr
  end select

  !-----------------------------------------------------------
  ! Define a maximal number of iteration
  maxcycle = 10000
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! Define the maximal number of generator / cells (for memory allocation purposes)
  print*,' Max number of generators : '
  max_number_generators = 50000 ;  print*,max_number_generators
  !----------------------------------------------------------

  !-----------------------------------------------------------
  !
  !   G L O B A L    A L L O C A T I O N S
  !
  max_points = MAXP
  allocate( XYp(1:max_points,1:2), XYp0(1:max_points,1:2), WW_c(1:max_points) )
  allocate( Vp(1:100), Cp(1:max_points,1:2), Npa(1:max_points) )
  allocate( W_p(1:max_points), W_c(1:max_points) )

  print*,' --> Allocations performed'
  !----------------------------------------------------------

  !--------------------------------------------------------------------------------
  !
  !   I N I T I A L I S A T I O N     O F    P R O B L E M
  ! Init problem type
  call init_problem(  XYp, Vp, Cp, Npa, problem, npart0, xmin,xmax,ymin,ymax )
  npart = npart0
  ! Copy the initial generators
  XYp0 = XYp
  print*,' --> Initialisation performed npart=',npart0
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !  O P E N      S C R I P T   F I L E S
  open(103,file='script.gnu')
  write(103,*) 'reset'
  ! print*, "ATTTEEENTION - 1"
  write(103,*) 'xmax=',xmax
  ! print*, "ATTTEEENTION - 2"

  write(103,*) 'xmin=',xmin
  write(103,*) 'ymax=',ymax
  write(103,*) 'ymin=',ymin
  write(103,*) 'dx=',(xmax-xmin)/10.0_d
  write(103,*) 'dy=',(ymax-ymin)/10.0_d
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !     M E S H E R    O F    V O R O N O I      K I N D
  !
  ! 1- Deallocate the mesh before
  if( debug == 1 ) print*,'   -- > Dealloc all in Mesh '
  call dealloc_ini_all(Mesh)
  !
  ! 2- Copy generators into x, y arrays
  if( debug == 1 ) print*,'   -- > Copy generators  npart=',npart
  x(1:npart) = XYp(1:npart,1)
  y(1:npart) = XYp(1:npart,2)
  n  = npart
  ! ng : nb of generator on boundary
  if( ng < 0 ) then
    ng = 4
  end if
  !
  ! 3- Alloc arrays
  if( debug == 1 ) print*,'   -- > Alloc Ini   ng,n=',ng,n
  call alloc_ini(Mesh)



  icycle = 1
  open (12,file="energy")
  do while ( icycle < 400 )
    !
    ! 4- Create Voronoi mesh
    if( debug == 1 ) print*,'   -- > Make Voronoi   ng,n=',ng,n
    call make_voronoi
    !
    ! 4- Create Voronoi mesh into data structure Mesh%
    critical_length = 0.0_d  ;  critical_angle  = 0.0_d
    if( debug == 1 ) print*,'   -- > Voronoi to staggered'
    call voronoi_to_staggered(critical_length,critical_angle)

    XYp(1:Mesh%nc,1) = x(1:Mesh%nc)
    XYp(1:Mesh%nc,2) = y(1:Mesh%nc)

    ! x(1:npart) = XYp(1:npart,1)
    ! y(1:npart) = XYp(1:npart,2)

    allocate( Mesh%X_n_n(1:Mesh%nn), Mesh%Y_n_n(1:Mesh%nn) )
    allocate( Mesh%X_c(1:Mesh%nc), Mesh%Y_c(1:Mesh%nc) )

    !
    !
    !      E N D    O F    M E S H E R
    !-----------------------------------------------------


    !------------------------------------------------------
    ! Compute subcell number
    Mesh%ncn = 0
    do i = 1,Mesh%nc
      Mesh%ncn = Mesh%ncn + Mesh%c_l(i)
    end do
    !-------------------------------------------------------
    ! Allocation of (node) structure for node positions
    !-------------------------------------------------------

    !------------------------------------------------------------------
    ! Copy node position
    Mesh%X_n_n(1:Mesh%nn) = xn(1:Mesh%nn)
    Mesh%Y_n_n(1:Mesh%nn) = yn(1:Mesh%nn)
    !-------------------------------------------------------------------

    !-------------------------------------------------------
    ! Allocation of (centroid) arrays for centroid position
    ! Compute centroids
    !call compute_centroids( Mesh )
    call compute_centroids_via_triangles( Mesh, XYp )

    if (icycle == 1 ) then
      energy0=0
        do p=5,Mesh%nc
          energy0=energy0+sqrt((x(p)-Mesh%X_c(p))**2 + (y(p)-Mesh%Y_c(p))**2 )
        enddo
    endif
    !-------------------------------------------------------

    !--------------------------------------------------------------------
    !
    !  C R E A T I O N       O F      F I L E S
    !

    write(cnum,*) icycle + maxcycle

    cnum     = adjustl(cnum)
    ! OUTPUT file part      ==> File of generators
    call create_part_file(filename,icycle,stringg,cnum,npart,XYp)
    print*,'   Generator file into ',filename

    ! OUTPUT FILE mesh ==> File mesh.#   polygons
    call create_mesh_file(filename2,icycle,stringg,cnum,Mesh,XYp,WW_c)
    print*,'   Mesh file into ',filename2
    ! OUTPUT FILE mesh ==> File tria.#   triangles
    ! call create_trimesh_file(filename2,icycle,stringg,cnum,Mesh,XYp,WW_c)
    ! print*,'   TriMesh file into ',filename2

    ! OUTPUT FILE mesh ==> File centr.#  centroids
    filename3 = 'centr'
    call create_centroid_file(filename3,icycle,stringg,cnum,Mesh )
    print*,'   Centroid file into ',filename3

    ! SCRIPT FOR GNUPLOT
    write(103,*) '  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy+dy] "',&
    & filename,'" t "Generator#',Mesh%nc,'" w p pt 5 ps 1 lc 3,"', &
    & filename2,'" t "Mesh" w l lt 1 lc 3,"', &
    & filename3,'" t "Centroids#',Mesh%nc,'" w p pt 6 ps 1 lc 2,"', &
    & filename3,'" t "time#',icycle,'" w p pt 6 ps 1 lc 2'

    write(103,*) ' pause 0.025'

    !------------------------------------------------------------------

    !----------------------------------------
    ! Elephant place

    !-------------------------------------------------------------------

    !--------------------------------------------------------------------
    !
    !  C R E A T I O N       O F      F I L E S
    !
    !!icycle = 1
    !if( in==inn ) then
    ! icycle = in
    ! write(cnum,*) icycle + maxcycle
    ! cnum     = adjustl(cnum)
    ! ! OUTPUT file part      ==> File of generators
    ! call create_part_file(filename,icycle,stringg,cnum,npart,XYp)
    ! print*,'   Generator file into ',filename
    !
    ! ! OUTPUT FILE mesh ==> File mesh.#
    ! call create_mesh_file(filename2,icycle,stringg,cnum,Mesh,XYp,WW_c)
    ! print*,'   Mesh file into ',filename2
    ! ! SCRIPT FOR GNUPLOT
    ! write(103,*) '  p [xmin-dx:xmax+dx][ymin-dy:ymax+dy+dy] "',&
    ! & filename,'" t "Generator#',Mesh%nc,'" w p pt 5 ps 1 lc 3,"', filename2,'" t "Mesh',in,'," w l lt 1 lc 3 '
    !------------------------------------------------------------------
    icycle = icycle + 1
    energy = 0.0
    do p = 5,Mesh%nc
      energy=energy+sqrt((x(p)-Mesh%X_c(p))**2 + (y(p)-Mesh%Y_c(p))**2)
      x(p) = Mesh%X_c(p)
      y(p) = Mesh%Y_c(p)
      ! call random_number(XYp(p,1))
      ! call random_number(XYp(p,2))
      ! XYp(p, 1) = Mesh%X_c(p)
      ! XYp(p, 2) = Mesh%Y_c(p)
    end do

    print*, icycle,energy,energy0,(energy/energy0)*100.0_d,"%"
    write(12,*) icycle,energy
    if (energy/energy0*100.0_d<1) then
      exit
    endif
    call odour_transmission(Mesh,icycle, W_c, W_p)




   !  if ( modulo(icycle, 20) == 0) then
   !  ! if (energy/energy0*100.0_d<5) then
   ! !
   !  npart = Mesh%nc+5
   !  do i = 1, 5
   !    call random_number(xx)
   !    call random_number(yy)
   !    x(mesh%nc+i) = xx/8_d + 0.4325_d
   !    y(mesh%nc+i) = yy/8_d + 0.4325_d
   !  end do
   ! end if
    n = npart
    ng = 4
    !-------------------------------------------------------
    ! Deallocation of  (node) structure for node positions
    !-------------------------------------------------------
    deallocate(Mesh%X_n_n, Mesh%Y_n_n)
    deallocate(Mesh%X_c, Mesh%Y_c)
  end do
  !----------------------------------------
  !   C L O S I N G    A R E A
  ! Close script files
  close(103)
  close(12)



!call potential(t,t0,t1,a,b,pot)
  !open(13,file="potentiel")
  !do p = 0,30

    ! pot = potential(p,10,20,0.5_d,1.0_d)

    ! write(13,*) p,pot
  ! end do

  ! close(13)
  !----------------------------------------

  !------------------------------------------------------------
  !
  !  F I N A L    D E A L L O C A T I O N
  !

  deallocate( XYp, XYp0, Vp, Cp, Npa, WW_c, W_p, W_c )
  !------------------------------------------------------------

  print*,'  OPEN GNUPLOT AND type >  load "script.gnu" '

end program VORONOI
