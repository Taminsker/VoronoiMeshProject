module bibli_init

  use mesh_data

  implicit none

contains

  ! INITIALISATION OF PROBLEM TYPE
  ! 0: random position of generators
  ! 1: uniform
  ! 2: -
  ! 3: -
  ! 4: Nothing
  ! 5: Annulus shape
  ! 6: full disk
  ! 7: -
  ! default: read a file of generators or mesh

  subroutine init_problem( XYp, Vp, Cp, Np, problem, nnature, xmin,xmax,ymin,ymax )
    implicit none
    ! Initialize the problem to be run
    !-------------
    ! Parameters
    real*8, dimension(:,:),  intent(inout) :: XYp
    real*8, dimension(:),    intent(inout) :: Vp
    integer,             dimension(:,:),  intent(inout) :: Cp
    integer,             dimension(:),    intent(inout) :: Np
    integer,                              intent(inout) :: problem,nnature
    real*8,                  intent(inout)    :: xmin,xmax,ymin,ymax
    !-------------
    character(30) :: filenamee
    integer :: i,j,ncells,npart,nx,ny, nppc, p,q, pp, ii, nn, nnp, sym,qq,mm,m
    real*8 :: dx, dy, drad,rad, ddx, ddy, xmn,xmx,ymn,ymx, xx,yy,pi,R,dR,R0,R1
    real*8 :: theta1,theta0,theta, dtheta
    !
    print*, '   ----> Enter Initialization problem #' , problem
    xmn = xmin; xmx = xmax
    ymn = ymin; ymx = ymax
    !
    ! Init the Generators
    !
    i = 1  ;   j = 1
    select case( problem )
    case(0)
      ! Random
      p = 0
      nppc = 2500
      print*,'   How many generators ?'
      read*,nppc
      nppc = nppc + 4
      print*,' Add 4 generators at the four box corners'
      ng = 4
      XYp(1,1) = xmin; XYp(1,2) = ymin
      XYp(2,1) = xmax; XYp(2,2) = ymin
      XYp(3,1) = xmax; XYp(3,2) = ymax
      XYp(4,1) = xmin; XYp(4,2) = ymax
      p = 4
      qq = p+1
      do q = qq,nppc
        p = p + 1
        call random_number(XYp(p,1))
        XYp(p,1) = xmn + abs(xmx-xmn) * XYp(p,1)
        call random_number(XYp(p,2))
        XYp(p,2) = ymn + abs(ymx-ymn) * XYp(p,2)
        Np(p) = 1
      end do
      !print*,' Add several generator around a specific point (1:YES, /=1 NO)'
      !read*,i
      !if( i==1 ) then
      !   print*,'  ==> Point to refine (x,y) and number of extra generator' ; read*,xx,yy,ng
      !   qq = p + 1
      !   nppc= nppc +  ng
      !   xmx = 0.1d0; ymx = 0.1d0
      !   do q = qq,nppc
      !      p = p + 1
      !      call random_number(XYp(p,1))
      !      XYp(p,1) = xx + abs(xmx-xmn) * XYp(p,1)
      !      call random_number(XYp(p,2))
      !      XYp(p,2) = yy + abs(ymx-ymn) * XYp(p,2)
      !      Np(p) = 1
      !   end do
      !end if
      npart = p
      ng = -1
      print*,'   ..Generators OK p=',p
    case(1)
      ! Uniform
      p = 0
      nppc = 250
      print*,'   How many generators (ex: 125,2500)?'
      read*,nppc
      dx = (xmx-xmn) / (sqrt(real(nppc)) )
      dy = (ymx-ymn) / (sqrt(real(nppc)) )
      ii = sqrt(real(nppc))
      print*,' Delta x, Delta y=',dx,dy,' Nc x Nc=',ii,'x',ii
      pp = 0
      nppc = nppc + 4
      print*,' Add 4 generators at the four box corners'
      ng = 4
      XYp(1,1) = xmin; XYp(1,2) = ymin
      XYp(2,1) = xmax; XYp(2,2) = ymin
      XYp(3,1) = xmax; XYp(3,2) = ymax
      XYp(4,1) = xmin; XYp(4,2) = ymax
      pp = 4
      do q = 1,ii
        do p = 1,ii
          pp = pp + 1
          if( mod(q,2)==0 ) then
            XYp(pp,1) = xmn + dx/2.d0 + (p-1)*dx
            XYp(pp,2) = ymn + dy/2.d0 + (q-1)*dy
          else
            XYp(pp,1) = xmn + dx/2.d0 + (p-1)*dx -dx/5.d0
            XYp(pp,2) = ymn + dy/2.d0 + (q-1)*dy -dy/5.d0
          end if
          Np(pp) = 1
        end do
      end do
      p = pp
      npart = p
      ng = -1
      print*,'   ..Generators OK p=',p
    case(2)
      !
    case(3)
      !
    case(4)
      !
    case(5)
      !
    case(6)
      !
    case(7)
      !
    case default
      ! Read a file of generators
      filenamee = 'T.mesh'
      open(11,file=filenamee)
      read(11,*) npart, ng
      do i = 1,npart
        read(11,*) XYp(i,1),XYp(i,2)
      end do
      close(11)
      print*,'   ..Generators OK p=',npart
    end select

    ! Final computational bounds
    xmin = minval(XYp(:,1)); ymin = minval(XYp(:,2));
    xmax = maxval(XYp(:,1)); ymax = maxval(XYp(:,2));
    ! Final number of particles
    nnature = p
    !
    open(11,file='end_init.plt')
    do i = 1,npart
      write(11,*) XYp(i,1),XYp(i,2)
    end do
    close(11)
  end subroutine init_problem


  subroutine create_wcentroid_file(filename,string,cnum,Mesh )
    implicit none
    !-------------
    ! Parameters
    type(Mesh_struct),     intent(in) :: Mesh
    character(11),intent(inout)       :: string
    character(30) ,intent(inout)      :: cnum
    character(10),intent(inout) :: filename
    integer :: p
    !------------
    string  = 'wcen.'
    filename = trim(string)//trim(cnum)
    open(21,file=filename)
    do p = 1,Mesh%nc
      write(21,*) Mesh%X_c(p),Mesh%Y_c(p)
    end do
    close(21)
  end subroutine create_wcentroid_file

  subroutine create_centroid_file(filename3,icycle,stringg,cnum,Mesh )
    implicit none
    !-------------
    ! Parameters
    type(Mesh_struct),     intent(in) :: Mesh
    character(11),intent(inout)       :: stringg
    character(30) ,intent(inout)      :: cnum
    integer,intent(inout)                :: icycle
    character(10),intent(inout) :: filename3
    integer :: p
    !------------
    stringg  = 'cent.'
    filename3 = trim(stringg)//trim(cnum)
    open(21,file=filename3)
    do p = 1,Mesh%nc
      write(21,*) Mesh%X_c(p),Mesh%Y_c(p)
    end do
    close(21)
  end subroutine create_centroid_file

  subroutine create_part_file( filename, icycle,stringg,cnum,npart,XYp )
    implicit none
    !-------------
    ! Parameters
    real*8, dimension(:,:),  intent(in) :: XYp
    character(11),intent(inout)       :: stringg
    character(30) ,intent(inout)      :: cnum
    integer,intent(inout)                :: npart,icycle
    character(10),intent(inout) :: filename
    integer :: p
    !------------
    stringg  = 'part.'
    filename = trim(stringg)//trim(cnum)
    !print*,' Cycle =',icycle,' file=',filename
    open(21,file=filename)
    do p = 1,npart
      write(21,*) XYp(p,1),XYp(p,2)
    end do
    close(21)
  end subroutine create_part_file


  subroutine create_mesh_file( filename2,icycle,stringg,cnum,Mesh,XYp,Res)
    implicit none
    !-------------
    ! Parameters
    type(Mesh_struct),     intent(in) :: Mesh
    real*8, dimension(:,:),  intent(in) :: XYp
    real*8, dimension(:),  intent(in) :: Res
    character(11),intent(inout)       :: stringg
    character(30) ,intent(inout)      :: cnum
    integer,intent(inout)                :: icycle
    character(10),intent(inout) :: filename2
    integer :: i,j
    real*8 :: xx,yy
    !------------
    stringg  = 'mesh.'
    filename2 = trim(stringg)//trim(cnum)
    filename2 = trim(filename2)
    open(24,file=filename2)
    xx = 0.9d0
    xx = 1.0d0
    yy = 1.0d0 - xx
    do j = 1,Mesh%nc

      do i = 1,Mesh%c_l(j)
        write(24,*) xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1) ,&
        &        xx* yn( Mesh%cell_list(j,i) ) + yy * XYp(j,2) , Res(j)
      end do
      i=1
      write(24,*) xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1) ,&
      &        xx* yn( Mesh%cell_list(j,i) ) + yy * XYp(j,2) , Res(j)
      write(24,*) ''
      write(24,*) ''
    end do
    close(24)
  end subroutine create_mesh_file




end module bibli_init
