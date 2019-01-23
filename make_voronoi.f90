subroutine make_voronoi
  ! control subroutine for generation of Voronoi (Dirichlet) mesh
  ! using routines supplied by Misha. The routines have been altered
  ! so that they compile using an F90 compiler
  use mesh_data
  !use tridata
  !
  implicit none
  !
  integer :: rc,i
  integer, dimension(10) :: inf
  !
  integer, external :: iDirGen 
  !
  !
  !write(*,*) 'Generating Voronoi cells'
  !print*,'   ** Subroutine make_voronoi - n=',n
  !
  ! Now generate Voronoi (Dirichlet) mesh from points above
  !
  ! Build the triangulation
  !
  call Delon(rc,inf)
  if(rc.ne.0) then
     write(*,*) ' Error while triangulation. pt too close rc,inf(1),inf(2)=',rc, inf(1), inf(2),n
     write(*,*) ' X1,Y1=',x(inf(1)),y(inf(1))
     write(*,*) ' X2,Y2=',x(inf(2)),y(inf(2))
     write(*,*) ' Dist=',sqrt( (x(inf(1))-x(inf(2)))**2+(y(inf(2))-y(inf(1)))**2 )
     open(14)
     do i = 1,n
        write(14,*) x(i),y(i),i
     end do
     close(14)
     stop
  endif
  !
  ! Build boundary:
  !
  !print*,'Bound'
  call Bound
  !print*,'End Bound'
  !
  ! Build Dirichlet tesselation:
  !
  !print*,'iDirGen'
  rc=iDirGen()
  if(rc.ne.0) then
     write(*,*) ' rc=',rc
     stop
  endif
  !
  ! I required call gnuplot output routine as test
  !call grafD2
  !call grafD
  !stop
  !
  !print*,'   ** End Subroutine make_voronoi'
  !
end subroutine make_voronoi
