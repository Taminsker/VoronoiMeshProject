module bibli_init

  use mesh_data
  use our_module

  implicit none

contains


  subroutine compute_centroids_via_triangles( Mesh, XYp )
    implicit none
    type(Mesh_struct),      intent(inout) :: Mesh
    real*8, dimension(:,:), intent(inout) :: XYp
    integer :: ic,j,in,in1
    real*8  :: SurfC,SurfT,Xc,Yc,x1,y1,x2,y2,x0,y0,Xt,Yt
    do ic = 1, Mesh%nc
       Xc = 0.0d0
       Yc = 0.0d0
       SurfC = 0.0d0
       ! Loop over the nodes of cell i
       do j = 1,Mesh%c_l(ic)-1
          ! jth and (j+1)th nodes of cell number ic
          in  = Mesh%cell_list(ic,j)
          in1 = Mesh%cell_list(ic,j+1)
          ! Edge [j,j+1], coordinates (x1,y1), (x2,y2)
          x1 = Mesh%X_n_n(in)    ; y1 = Mesh%Y_n_n(in)
          x2 = Mesh%X_n_n(in1)   ; y2 = Mesh%Y_n_n(in1)
          ! Generator position
          x0 = XYp(ic,1) ; y0 = XYp(ic,2)
          ! Triangle T=(X0,X1,X2)
          SurfT =  abs(0.5d0*( (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0) ))
          ! Centroid of triangle T : Xt,Yt
          Xt = 1.0d0/6.0d0*( (y1-y0)*(x1**2+x1*x0+x0**2) + &
               &             (y2-y1)*(x2**2+x2*x1+x1**2) + &
               &             (y0-y2)*(x0**2+x0*x2+x2**2) )
          Yt =-1.0d0/6.0d0*( (x1-x0)*(y1**2+y1*y0+y0**2) + &
               &             (x2-x1)*(y2**2+y2*y1+y1**2) + &
               &             (x0-x2)*(y0**2+y0*y2+y2**2) )
          ! The formula above are equivalent to this one !!!!
          !Xt = 1.0d0/3.0d0*( x0+x1+x2 )
          !Yt = 1.0d0/3.0d0*( y0+y1+y2 )
          write(101,*) Xt,Yt,SurfT
          ! Cell centroid : cumulate triangle centroids
          Xc = Xc + Xt
          Yc = Yc + Yt
          ! Cell surface by cumulative triangle surfaces
          SurfC = SurfC + SurfT
       end do
       j = Mesh%c_l(ic)
       in  = Mesh%cell_list(ic,j)
       in1 = Mesh%cell_list(ic,1)
       ! Edge [j,j+1], coordinates (x1,y1), (x2,y2)
       x1 = Mesh%X_n_n(in)    ; y1 = Mesh%Y_n_n(in)
       x2 = Mesh%X_n_n(in1)   ; y2 = Mesh%Y_n_n(in1)
       ! Triangle T=(X0,X1,X2)
       SurfT =  abs(0.5d0*( (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0) ))
       ! Centroid of triangle T : Xt,Yt
       Xt = 1.0d0/6.0d0*( (y1-y0)*(x1**2+x1*x0+x0**2) + &
            &             (y2-y1)*(x2**2+x2*x1+x1**2) + &
            &             (y0-y2)*(x0**2+x0*x2+x2**2) )
       Yt =-1.0d0/6.0d0*( (x1-x0)*(y1**2+y1*y0+y0**2) + &
            &             (x2-x1)*(y2**2+y2*y1+y1**2) + &
            &             (x0-x2)*(y0**2+y0*y2+y2**2) )
       ! The formula above are equivalent to this one !!!!
       !Xt = 1.0d0/3.0d0*( x0+x1+x2 )
       !Yt = 1.0d0/3.0d0*( y0+y1+y2 )
       ! Cell centroid : cumulate triangle centroids
       Xc = Xc + Xt
       Yc = Yc + Yt
       ! Cell surface by cumulative triangle surfaces
       SurfC = SurfC + SurfT
       ! Scale by cell surface
       Xc = Xc / SurfC
       Yc = Yc / SurfC

       ! Store results
       Mesh%X_c(ic) = Xc
       Mesh%Y_c(ic) = Yc

    end do
  end subroutine compute_centroids_via_triangles


  subroutine compute_centroids( Mesh )
    implicit none
    type(Mesh_struct),intent(inout) :: Mesh
    integer :: i
    real*8  :: Sc
    do i = 1, Mesh%nc
       ! This call compute the cell surface
       Sc          = cell_center_massM(Mesh,i,0.0d0,0.0d0,1.0d0)
       Mesh%X_c(i) = cell_center_massM(Mesh,i,1.0d0,0.0d0,0.0d0) / Sc
       Mesh%Y_c(i) = cell_center_massM(Mesh,i,0.0d0,1.0d0,0.0d0) / Sc
    end do
  end subroutine compute_centroids

  function cell_center_massM(Mesh,ic,a,b,c)
    implicit none
    !
    ! f = a.x + b.y + c
    !
    ! (x1,y1),(x2,y2),(x3,y3),(x4,y4) is the counterclock wise oriented boundary
    ! of a corner
    type(Mesh_struct),intent(inout) :: Mesh
    real*8             :: cell_center_massM
    real*8,intent(in)  :: a,b,c
    integer,intent(in) :: ic
    !
    real*8 :: xx,yy,uu,vv,line_integral2,a6,b6,c2
    integer             :: i
    !
    line_integral2 = 0.0d0
    !
    a6 = a/6.0d0
    b6 = b/6.0d0
    c2 = c/2.0d0
    !
    !print*,' Centroid of cell',ic,a,b,c
    do i = 1,Mesh%c_l(ic)-1
       ! Find the points for performing the integral
       xx = Mesh%X_n_n( Mesh%cell_list(ic,i) )   ;  yy = Mesh%Y_n_n( Mesh%cell_list(ic,i) )
       uu = Mesh%X_n_n( Mesh%cell_list(ic,i+1) ) ;  vv = Mesh%Y_n_n( Mesh%cell_list(ic,i+1) )
       !
       !print*,i,Mesh%c_l(ic), Mesh%cell_list(ic,i), Mesh%cell_list(ic,i+1)
       line_integral2 = line_integral2 &
            &  + a6 * (vv-yy)*(uu**2 + uu*xx + xx**2)  &
            &  - b6 * (uu-xx)*(vv**2 + vv*yy + yy**2)  &
            &  + c2 * (vv-yy)*(uu+xx)
    end do
    i = Mesh%c_l(ic)
    xx = Mesh%X_n_n( Mesh%cell_list(ic,i) )   ;  yy = Mesh%Y_n_n( Mesh%cell_list(ic,i) )
    uu = Mesh%X_n_n( Mesh%cell_list(ic,1) )   ;  vv = Mesh%Y_n_n( Mesh%cell_list(ic,1) )
    !print*,i,Mesh%c_l(ic), Mesh%cell_list(ic,i), Mesh%cell_list(ic,1)
    line_integral2 = line_integral2 &
         &  + a6 * (vv-yy)*(uu**2 + uu*xx + xx**2) &
         &  - b6 * (uu-xx)*(vv**2 + vv*yy + yy**2)  &
         &  + c2 * (vv-yy)*(uu+xx)
    !
    cell_center_massM = line_integral2
    !print*,cell_center_massM
    !
  end function cell_center_massM



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
    integer :: i,j,ncells,npart,nx,ny, nppc, p,q, pp, ii, nn, nnp, sym,qq,mm,m,ne
    real*8 :: dx, dy, drad,rad, ddx, ddy, xmn,xmx,ymn,ymx, xx,yy,pi,R,dR,R0,R1
    real*8 :: theta1,theta0,theta, dtheta
    real*8 :: s1x,s1y,s2x,s2y,s1l,s2l
    !
    print*,'   ----> Enter Initialization problem #',problem
    print*,' PROBLEM#',problem
    xmn = xmin; xmx = xmax
    ymn = ymin; ymx = ymax
    !
    ! Init the Generators
    !
    i = 1  ;   j = 1
    select case( problem )
    case(22)
       print*,' Elephant/generator initialisation, number of elephants?'
       read*,ne
       nppc = ne + 4 + 2*100
       ! Number of elephants
       ng = 4
       XYp(1,1) = xmin; XYp(1,2) = ymin
       XYp(2,1) = xmax; XYp(2,2) = ymin
       XYp(3,1) = xmax; XYp(3,2) = ymax
       XYp(4,1) = xmin; XYp(4,2) = ymax
       p = 4
       Np(1:p) = 1
       ! Generators on S1=(1/4,3/4), S2=(3/4,1/4)
       s1x = 0.25d0 ; s1y = 0.75d0 ; s1l = 0.05d0
       s2y = 0.25d0 ; s2x = 0.75d0 ; s2l = 0.075d0
       ! Source S1
       q = 0
       p = p + 1
       do
          call random_number(XYp(p,1))
          call random_number(XYp(p,2))
          if(    XYp(p,1)<= s1x+s1l .and. XYp(p,1)>= s1x-s1l .and. &
               & XYp(p,2)<= s1y+s1l .and. XYp(p,2)>= s1y-s1l ) then
             Np(p) = 1
             q = q + 1
             p = p + 1
          end if
          if( q==100 ) exit
       end do
       ! Source S2
       q = 0
       do
          call random_number(XYp(p,1))
          call random_number(XYp(p,2))
          if(    XYp(p,1)<= s2x+s2l .and. XYp(p,1)>= s2x-s2l .and. &
               & XYp(p,2)<= s2y+s2l .and. XYp(p,2)>= s2y-s2l ) then
             Np(p) = 1
             q = q + 1
             p = p + 1
          end if
          if( q==100 ) exit
       end do
       ! Generators everywhere
       p = p - 1
       qq = p
       do q = qq,nppc
          p = p + 1
          call random_number(XYp(p,1))
          call random_number(XYp(p,2))
          Np(p) = 2
       end do
       npart = p
       ng = -1
       print*,' Number of elephants =',ne
       print*,'   ..Generators OK p=',p
    case(0)
       ! Random
       p = 0
       nppc = 2500
       print*,'   How many elephants?'
       read*,nppc
       nppc = nppc + 4 + 2 *100
       print*,' Add 4 generators at the four box corners'
       ng = 4
       XYp(1,1) = xmin; XYp(1,2) = ymin
       XYp(2,1) = xmax; XYp(2,2) = ymin
       XYp(3,1) = xmax; XYp(3,2) = ymax
       XYp(4,1) = xmin; XYp(4,2) = ymax
       call zone(Mesh, XYp, 0.25_d, 0.25_d, 0.1_d, 0.75_d, 0.75_d, 0.1_d)
       ! p = 206 !4
       p = 4+2*100
       ! nppc = nppc + 200
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
       !filenamee = 'T.mesh'
       filenamee = 'world.txt'
       open(11,file=filenamee)
       read(11,*) npart, ng
       do i = 1,npart
          read(11,*) XYp(i,1),XYp(i,2)
       end do
       close(11)
       p = npart
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


subroutine create_wcentroid_file(filename3,icycle,stringg,cnum,Mesh )
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
    stringg  = 'wcen.'
    filename3 = trim(stringg)//trim(cnum)
    open(21,file=filename3)
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
    do p = 205,npart
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



  subroutine create_trimesh_file( filename2,icycle,stringg,cnum,Mesh,XYp,Res)
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
    integer :: i,j,in,ic1,ic2,ic3
    !------------
    stringg  = 'tria.'
    filename2 = trim(stringg)//trim(cnum)
    filename2 = trim(filename2)
    open(24,file=filename2)
    do in = 1,Mesh%nn
      !a utiliser pour le W_p et W_c
       if( Mesh%n_l(in)>=3 ) then
          ic1 = Mesh%node_list(in,1)
          ic2 = Mesh%node_list(in,2)
          ic3 = Mesh%node_list(in,3)
          write(24,*) XYp(ic1,1),XYp(ic1,2)
          write(24,*) XYp(ic2,1),XYp(ic2,2)
          write(24,*) XYp(ic3,1),XYp(ic3,2)
          write(24,*) XYp(ic1,1),XYp(ic1,2)
          write(24,*) ''
          write(24,*) ''
       end if
    end do
    close(24)
  end subroutine create_trimesh_file



end module bibli_init
