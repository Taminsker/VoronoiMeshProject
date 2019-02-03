!********************************************************************
!                             DELON                                 *
!********************************************************************
!   Delon - main subroutine for Delaunay triangulation of a set of  *
! points. This subroutine build a convex hull of the point set.     *
! rc - return code of the subroutine:                               *
!     rc=0 - no error found                                         *
!     rc=1 - two points from the set are close to each other.       *
!            inf(1) and inf (2) - the numbers of points.            *
!     rc=2 - you try to use more points then possible. Change       *
!            parameter MAXP in the file TRICONST.INC and recompile  *
!            files. inf(1) - current number of points, inf(2) -     *
!            maximal available number of points.                    *
!     rc=3 - used number of points is less then 3. It is impossible *
!            to build triangulation for this set of point.          *
!            inf(1) - current number of points.                     *
!     rc>200 - internal error. Please contact with authors and send *
!            them file with initial data (coordinates of points).   *
!********************************************************************
subroutine Delon(rc,inf)
!use tridata
use mesh_data
double precision xRand
common/cRand/ xRand
!
integer rc,inf(10),io,n10,ip
integer iRand,iFindTri               ! Functions
!
integer iSet(0:MAXTRI)
!-------------- Local variables --------------------------------
integer i,j,nt0
!
!print*,'      ** Subroutine delon.f90 - n=',n
io=0
!
if(n.lt.3) then
 rc=3
 inf(1)=n
 return
end if
!
if(n.ge.MAXP) then
 rc=2
 inf(1)=n
 inf(2)=MAXP-4
 return
end if
!
xRand=13000000001.
eps=1.e-12*sqrt((xmax-xmin)**2+(ymax-ymin)**2)
! 
nt0=1

!
call TriInit              ! Big triangle (2n+2)
!

do 100 i=1,n
 iSet(i)=-1              ! Clear inform array
100 continue
!
n10=n/10                  ! Part for random insert
!
do 200 iP=1,n10

   ! Add on Raph on Dec 2009
   do
      !i=iRand()
      !call rRand(r)
      call random_number(r)
      rn=n
      ir=rn*r
      if(ir.le.0) ir=1
      if(ir.ge.n) ir=n
      i=ir
      if( i<=size(iSet) ) exit
   end do
 

 !print*,' iP,i=',iP,i,' n10=',n10
 !print*,' iSet(i)=',iSet(i),' Size(iSet)=',size(iSet)

 if(iSet(i).eq.0) goto 200
 iSet(i)=0
 !
 j=iFindTri(i,rc)        ! Find trian. with the point i
 !print*,'j=',j
 if(rc.ne.0) return      ! Some error in GetTri
 !
 call CheckTri(i,j,rc,inf)   ! Check the distance between points
 !print*,'i,j,rc,inf=',i,j,rc,inf
 if(rc.ne.0) return      ! Points to close each other
 !
 call GetList(i,j,nt0)   ! Build a list of triang. to be corrected
 !print*,i,j,nt0
 if(rc.ne.0) return
 !
 call ChangeTri(i)       ! Correct triangles 
 !
200 continue
!

do 300 i=1,n
 !
 if(nt.ge.1299) then
  j=0
 endif
 !        write(88,'(a,i5)') '   Point',i
 if(i.ge.1295/2) then
  j=0
 endif
 !
 if(iSet(i).eq.0) goto 300
 iSet(i)=0
 !
 j=iFindTri(i,rc)            ! Find trian. with the point i
 if(rc.ne.0) return          ! Some error in GetTri
 !
 call CheckTri(i,j,rc,inf)   ! Check the distance between points
 if(rc.ne.0) return          ! Points to close each other
 !
 call GetList(i,j,nt0)       ! Build a list of triang. to be correct
 if(rc.ne.0) return
 !
 call ChangeTri(i)           ! Correct triangles
300 continue
!
call GetFirst                ! Build references to first triangles
!
rc=0
!
!print*,'      ** End Subroutine delon.f90'
return
end
!
!********************************************************************
!                            CHECKTRI                               *
!********************************************************************
!   This procedure check if new point lies in the same position as  *
! some oter point. New point has the number I and triangle with     *
! this point has the number IT. We test three points that are in    *
! the vertixes of the triangle.                                     *
!********************************************************************
subroutine CheckTri(i,it,rc,inf)
use mesh_data 
!
integer i,it,rc,inf(10)
integer j,k
double precision xi,yi,xk,yk
!
xi=x(i)
yi=y(i)
!print*,' i=',i,' xi,yi=',xi,yi
!
do j=0,2                  ! For all vertexes of triangle IT
   k=tri(it)%v(j)
 xk=x(k)
 yk=y(k)
 !print*,'k=',k,' xk,yk=',xk,yk
 if((xk-xi)**2+(yk-yi)**2.le.eps) then
    print*,'k=',k,' xk,yk=',xk,yk,' xi,yi=',xi,yi,' d=',(xk-xi)**2+(yk-yi)**2,'eps=',eps
  inf(1)=i
  inf(2)=k
  rc=1
  return
 end if
end do
!
return
end
!********************************************************************
!   Initialization procedure for triangle generation.               *
!   Here we'll build very big triangle. All points must lies in it. *
!********************************************************************
subroutine TriInit
use mesh_data 
!-------------- Local variables --------------------------------
!
double precision dx,dy,dx2,dy2
integer ia1,ia2,ia3,ia4
!
RIGHT(0)=1
RIGHT(1)=2
RIGHT(2)=0
!
LEFT(0)=2
LEFT(1)=0
LEFT(2)=1
!
dx=999.*(xmax-xmin)
dy=999.*(ymax-ymin)
!
nt=1
!
tri(0)%v(2)=n+1           ! Inner big triangle (vertexes)
tri(0)%v(1)=n+2
tri(0)%v(0)=n+3
!
tri(0)%t(0)=1             ! Inner big triangle (neighbours)
tri(0)%t(1)=1
tri(0)%t(2)=1
!
tri(1)%v(0)=n+1           ! Outer big triangle (vertexes)
tri(1)%v(1)=n+2
tri(1)%v(2)=n+3
!
tri(1)%t(0)=0             ! Outer big triangle (neighbours)
tri(1)%t(1)=0
tri(1)%t(2)=0
!
dx2=dx+dx
dy2=dy+dy
!
x(n+1)=xmin-dx2           ! Coordinates of
y(n+1)=ymin-dy            !   3 fictiv points
x(n+2)=xmax+dx2
y(n+2)=ymin-dy
x(n+3)=xmin
y(n+3)=ymax+dy2
!
ia1=n+1
ia2=n+2
ia3=n+3
ia4=1
call GetV(ia1,ia2,ia3,ia4)
!
return
end
!********************************************************************
!  GetV  -  calculate the centre of triangle                        *
!        i1,i2,i3 - the numbers of the vertexes of the triangle     *
!        iv - the number of triangle                                *
!********************************************************************
subroutine GetV(i1,i2,i3,iv)
use mesh_data 
integer i1,i2,i3,iv

double precision xi1,yi1,xj1,yj1,xk1,yk1
double precision a0,a1,zn

xi1=x(i1)           ! Coordinates of the vertexes
yi1=y(i1)
xj1=x(i2)
yj1=y(i2)
xk1=x(i3)
yk1=y(i3)

a0=xj1*xj1-xi1*xi1+yj1*yj1-yi1*yi1
a1=xk1*xk1-xi1*xi1+yk1*yk1-yi1*yi1
zn=1./(2.*((yk1-yi1)*(xj1-xi1)-(yj1-yi1)*(xk1-xi1)))

xv(iv)=((yk1-yi1)*a0-(yj1-yi1)*a1)*zn
yv(iv)=((xi1-xk1)*a0-(xi1-xj1)*a1)*zn

return
end

!********************************************************************
!       Generator of pseudo-random number in the range 1 - N        *
!********************************************************************

integer function iRand()
use mesh_data 
double precision r,rn
integer ir
!
call rRand(r)
rn=n
ir=rn*r
if(ir.le.0) ir=1
if(ir.ge.n) ir=n
iRand=ir
!
return
end
!
!********************************************************************
!
subroutine rRand(r)
double precision dm37
double precision r
double precision x
common/cRand/ x
!
dm37=137438953472.
x=x*3125.
x=x-idint(x/dm37)*dm37
r=x/dm37
!
return
end
!
!********************************************************************
!                       FINDTRI          GETTRI                     *
!********************************************************************
!   These two functions find a triangle with the point number I.    *
!  rc - return code.                                                *
!********************************************************************
!
integer function iFindTri(i0,rc)
use mesh_data 
!
integer i0,rc
integer iGetTri
!
rc=0
!
if(nt.eq.1) then                          ! Only one triangle
  iFindTri=1
  return
endif
!
iFindTri=iGetTri(x(i0),y(i0),rc)
!
return
end
!
!********************************************************************
!
integer function iGetTri(x0,y0,rc)
use mesh_data 
!
double precision x0,y0
integer rc
!
integer jt,jtt,ii,i,j,ir,il,is
double precision x1,y1,x2,y2,x3,y3,Stri
logical FirstHere
!
FirstHere=.TRUE.
jt=nt
jtt=jt
!
do i=0,2                                  ! Find real point ..
 j=tri(jt)%v(i)                          ! .. in the triangle nt
 if(j.le.n) goto 1000
end do
1000 continue
!
x1=x(j)
y1=y(j)
!
ii=0
!
100 continue
!
i=ii
!
110 continue
!
if(tri(jt)%v(i).eq.j) then
 !
 101     ir=RIGHT(i)
 il=LEFT(i)
 !
 if(FirstHere) then
  is=tri(jt)%v(ir)                    ! Right neighbour of point i
  x3=x(is)                            ! Coordinares of this ..
  y3=y(is)                            !     .. neighbour
 endif
 !
 is=tri(jt)%v(il)                ! Left neighbour
 x2=x3                                 ! Right neighbour coordinates
 y2=y3
 x3=x(is)                              ! Left neighbour coordinates
 y3=y(is)
 !
 STri=(x0-x1)*(y0-y2)-(y0-y1)*(x0-x2)
 if(STri.lt.0) goto 201                ! To next triang. over Pj
 !
 STri=(x0-x3)*(y0-y1)-(y0-y3)*(x0-x1)  ! To next triang. over Pj
 if(STri.lt.0) goto 201
 !
 STri=(x0-x2)*(y0-y3)-(y0-y2)*(x0-x3)
 if(STri.ge.0) then                    ! We have found triangle!
  iGetTri=jt
  return
 end if
 !
 jt=tri(jt)%t(i)
 jtt=jt
 !
 do i=0,2
  if(tri(jt)%v(i).eq.is) then
   ii=LEFT(i)
   j=tri(jt)%v(ii)             ! New point Pj
   x1=x(j)
   y1=y(j)
   goto 100
  end if
 end do
 !
 rc=200
 iGetTri=0
 return
 !
endif
!
i=i+1
if(i.lt.3) goto 110
!
201 FirstHere=.FALSE.                         ! Don't use right triangle
jt=tri(jt)%t(ir)                          ! Next triangle
ii=0
if(jt.ne.jtt) goto 100
!
rc=201                                   ! Internal error
iGetTri=0                                ! Triangle not found
return
end
!
subroutine GetFirst
use mesh_data !tridata
!
integer i,j,it,ii
!
do j=0,n+3
 ftri(j)=0
end do
!
do it=1,nt
 do i=0,2
  j=tri(it)%v(i)
  !
  ii=LEFT(i)
  !
  if(ftri(j).eq.0 .or. tri(it)%t(ii).lt.1) ftri(j)=it
 end do
end do
!
return
end
!
!-------------------------------------------------------------------
!
subroutine GetList(i0,t0,nt0)
use mesh_data 
!
integer i0,t0,nt0
!
integer v,tt,i1,i2,ii,it
double precision dx,dy,ri,r,x0,y0
!
x0=x(i0)
y0=y(i0)
!
l1(0)=t0
l2(0)=t0
n1=1
n2=1
!
i1=0
!
do while(i1.lt.n1)                ! For all elements in array L1
 tt=l1(i1)                       !
 do ii=0,2
  it=tri(tt)%t(ii)
  if(it.eq.0) goto 1
  if(it.lt.nt0) goto 1    ! Pass fictitious triangle
  do i2=0,n2-1            ! Find IT in the list L2
   if(l2(i2).eq.it) goto 1     ! Goto 1 if it is in the L2
  end do
  !
  l2(n2)=it               ! Put IT in the list L2
  n2=n2+1
  !
  v=tri(it)%v(0)                ! One of vertexes of triangle IT
  !
  dx=xv(it)-x(v)
  dy=yv(it)-y(v)
  r=dx**2+dy**2                 ! Radius of triangle IT
  !
  dx=xv(it)-x0
  dy=yv(it)-y0
  ri=dx**2+dy**2                ! Distance to point I0
  !
  if(ri.lt.r) then        ! Put IT in the list L1
   l1(n1)=it
   n1=n1+1
  endif
  !
  1   continue
 end do
 i1=i1+1
end do
!
return
end
!
!-------------------------------------------------------------------------
!
subroutine ChangeTri(i0)
use mesh_data 
integer i0
!
integer i,is,i1,i2,j1,j2,j3,iv,iv0
integer it,tt,pt,ii
!
if(nt.eq.1) then                  ! Build first triangle
 l2(0)=n+1
 l2(1)=n+2
 l2(2)=n+3
 n2=3
 !
 l3(0)=0
 l3(1)=0
 l3(2)=0
 !
 goto 2
endif
!
! Find a triangle in L1 a neibhour for which is absent in L1
!
do i1=n1-1,0,-1
 tt=l1(i1)                       ! Triangle from L1
 do i=0,2                        ! For all neighbours
  i2=n1
  it=tri(tt)%t(i)               ! IT - neighbour for TT
  if(it.ne.0) then
   do i2=0,n1-1                ! Find triangle IT ...
    if(l1(i2).eq.it) goto 10  !   ... in the list L1
   end do
  end if                  ! We didn't find
  !
  is=RIGHT(i)              ! Right neighbour index
  !
  iv0=tri(tt)%v(is)       ! Stop point
  pt=tt
  !
  is=LEFT(i)                    ! Left neighbour index
  !
  iv=tri(tt)%v(is)
  l2(0)=iv
  l3(0)=it
  n2=1
  goto 1
  !
  10     continue
 end do
end do
!
write(*,*)' Internal error 01'
stop                              ! We didn't find triangle !!! ???
!
1 continue
!
! Go around vertex IV and find triangle that is absent in L1. Start from PT
!
ii=0
!
3 do i=ii,2
 if(tri(pt)%v(i).eq.iv) then     ! We have find vertex IV in PT
  i2=n1
  !
  ii=LEFT(i)
  !
  it=tri(pt)%t(ii)              ! IT - the neighbour of triangle TT
  if(it.ne.0) then
   do i2=0,n1-1                ! Find triangle IT in list L1
    if(l1(i2).eq.it) then
     pt=it
     goto 13
    endif
   end do
  endif
  !                                ! Didn't find
  ii=RIGHT(i)
  !
  iv=tri(pt)%v(ii)        ! New vertex of polyhedra
  l2(n2)=iv
  l3(n2)=it                     ! Triangle across segment of boundary
  n2=n2+1
  if(iv.ne.iv0) goto 3
  !
  goto 13
 endif
end do
!
write(*,*)' Internal error 02'
stop                              ! Error. Didn't find IV in triangle PT
!
13 if(iv.ne.iv0) goto 1
!
2 l1(n1  )=nt+1
l1(n1+1)=nt+2
n1=n1+2
nt=nt+2
!
! Delete triangles from L1 and build new ones - with vertex in point I0
!
do i1=0,n1-1
 tt=l1(i1)
 !
 j1=i0                           ! Three vertexes of new triangle
 tri(tt)%v(0)=j1
 !
 ii=i1-1
 if(i1.eq.0) ii=n2-1
 j2=l2(ii)
 tri(tt)%v(1)=j2
 !
 j3=l2(i1)
 tri(tt)%v(2)=j3
 !
 it=l3(i1)                       ! Three neighbous of new triangle
 tri(tt)%t(0)=it
 !
 ii=i1+1
 if(ii.eq.n1) ii=0
 tri(tt)%t(1)=l1(ii)
 !
 ii=i1-1
 if(i1.eq.0) ii=n1-1
 tri(tt)%t(2)=l1(ii)
 !
 do i=0,2                        ! Find vertex J2 in the neighbouring
  if(tri(it)%v(i).eq.j2) then   !                             triangle
   !
   ii=RIGHT(i)
   !
   tri(it)%t(ii)=tt
   goto 14
  end if
 end do
 !
 write(*,*)' Internal error 03'
 stop                      ! Didn't found
 !
 14   call GetV(j1,j2,j3,tt)
 !
end do
!
return
end
