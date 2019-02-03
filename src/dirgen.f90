!**********************************************************************
!
integer function iAddNeib(iv)
use Mesh_data 
integer iv
dNeib(np)%neib=iv
dNeib(np)%next=np+1
np=np+1
iAddNeib=0
if(np.ge.MAXDIR) iAddNeib=2910
end
!
!**********************************************************************
!
integer function iDirGen()
use Mesh_data
integer ip,ip1,ip2,iv,it,it0
integer i,j,rc
!
integer iAddNeib,iGetVerN,iSgCorr,iVertDef    ! Functions
!
! print*,' IDIRGEN 0'
np=0
do ip=1,ng
 ip1=ip+1
 if(ip.EQ.ng) ip1=1
 ip2=ip-1
 if(ip2.EQ.0) ip2=ng
 it0=ftri(ip)
 it=it0
 dRef(ip)=np
 !
 if(iAddNeib(0).ne.0 .OR. iAddNeib(0).ne.0) then
  iDirGen=2910
  return
 end if
 !
 !  1. Find triangle with side (ip,ip1) (ip1 should be right for ip)
 !
 !print*,' IDIRGEN 1'
 it0=it
 10   continue
 i=iGetVerN(it,ip)
 j=RIGHT(i)
 if(tri(it)%v(j).EQ.ip1) goto 20
 it=tri(it)%t(j)
 if(it.ne.it0) goto 10
 !
 iDirGen=2913
 return
 !
 ! 2. Go around pointeger ip and write all connected points as neighbours
 !
 !print*,' IDIRGEN 2'
 20   it0=it
 30   continue
 i=iGetVerN(it,ip)
 j=RIGHT(i)
 iv=tri(it)%v(j)
 if(iAddNeib(iv).ne.0) then
  iDirGen=2910
  return
 end if
 if(iv.EQ.ip2) goto 40
 it=tri(it)%t(j)
 if(it.NE.it0) goto 30
 !
 iDirGen=2913
 return
 !
 40   dNeib(np-1)%next=dRef(ip)
 !
end do
!
!----------------------------------------------------------------------
!
! print*,' IDIRGEN 3'
do ip=ng+1,n
 dRef(ip)=np
 it0=ftri(ip)
 it=it0
 50   continue
 i=iGetVerN(it,ip)
 j=RIGHT(i)
 iv=tri(it)%v(j)
 if(iAddNeib(iv).ne.0) then
  iDirGen=2910
  return
 end if
 it=tri(it)%t(j)
 if(it.NE.it0) goto 50
 !
 dNeib(np-1)%next=dRef(ip)
 !
end do
!
! print*,' IDIRGEN 4'
do dPusto=np,MAXDIR
 dNeib(np)%next=np+1
end do
dNeib(MAXDIR)%next=0
!
rc=iSgCorr()
if(rc.NE.0) then
 iDirGen=rc
 return
end if
!
rc=iVertDef()
if(rc.NE.0) then
 iDirGen=rc
 return
end if
!
iDirGen=0
end
!
!**********************************************************************
!
integer function iVertDef()
use Mesh_data
integer i,j,jend,j0
integer p1,p2,pi,pk,i0,i1
double precision xx,yy,xi,yi,xj,yj,xk,yk,x1,y1,x2,y2,xs,ys,qq,a1,a2
!
!---------------------------------------------
!print*,size(dRef),dRef(:600)
!print*,n
do i=1,n
 xi=x(i)
 yi=y(i)
 j0=dRef(i)
 i1=dNeib(j0)%neib
 jend=dNeib(j0)%next
 j=jend
 !print*,i,':j0,i1,jend=',j0,i1,jend,' i0,i1=',i0,i1
 !
 10     i0=i1
 i1=dNeib(j)%neib
 !print*,i,':j0,i1,jend=',j0,i1,jend,' i0,i1=',i0,i1
 if(i0.EQ.0 .and. i1.EQ.0) then
  xx=xi
  yy=yi
  goto 30
 else if(i0.EQ.0) then
  p1=i
  p2=i+1
  if(i.EQ.ng) p2=1
  pi=i1
  pk=i
  goto 20
 else if(i1.EQ.0) then
  p1=i-1
  if(p1.eq.0) p1=ng
  p2=i
  pi=i0
  pk=i
  goto 20
 else
  xk=x(i0)
  yk=y(i0)
  xj=x(i1)
  yj=y(i1)
  qq=1./(2.*((yk-yi)*(xj-xi)-(yj-yi)*(xk-xi)))
  a1=xj*xj-xi*xi+yj*yj-yi*yi
  a2=xk*xk-xi*xi+yk*yk-yi*yi
  xx=((yk-yi)*a1-(yj-yi)*a2)*qq
  yy=((xi-xk)*a1-(xi-xj)*a2)*qq

!print*,xx,i0,i1,x(i0),y(i0),x(i1),y(i1)
!read*
  goto 30
 endif
 20     x1=x(p1)
 y1=y(p1)
 x2=x(p2)
 y2=y(p2)
 xj=x(pi)
 yj=y(pi)
 xk=x(pk)
 yk=y(pk)
 xs=0.5*(xj+xk)
 ys=0.5*(yj+yk)
 qq=(x2-x1)*(xk-xj)+(y2-y1)*(yk-yj)
 if(qq.EQ.0.) then
  iVertDef=2913
  return
 end if
 qq=((x1-xs)*(y2-ys)-(x2-xs)*(y1-ys))/qq
 xx=xs+(yk-yj)*qq
 yy=ys-(xk-xj)*qq
 !
 30     xVert(j0)=xx
 yVert(j0)=yy
 j0=j
 j=dNeib(j0)%next

!print*,xx,yy,j0,j

 if(j.NE.jend) goto 10
 !
end do
!
iVertDef=0
end
!
!**********************************************************************
!
subroutine DelNeib(i,prev,curr)
use Mesh_data
integer i,prev,curr
integer next
!
next=dNeib(curr)%next
dNeib(curr)%next=dPusto
dPusto=curr
dNeib(prev)%next=next
if(dRef(i).EQ.curr) dRef(i)=prev
end
!
!**********************************************************************
!
integer function iDelSide(it,ip)
use Mesh_data
integer it,ip
double precision x0,y0,x1,y1,dx,dy,d2
common/dDir2/ x0,y0,x1,y1,dx,dy,d2
integer i0,i1,i2
integer i,iend,iold,ii
integer neib0,neib1
double precision rq
!
i0=tri(it)%v(RIGHT(ip))
i1=tri(it)%v(LEFT(ip))
i2=tri(it)%v(ip)
!
iold=dRef(i0)
iend=dNeib(iold)%next
i=iend
!
10   if(dNeib(i)%neib.EQ.i1) then
 call DelNeib(i0,iold,i)
 goto 20
end if
iold=i
i=dNeib(iold)%next
if(i.ne.iend) goto 10
!
iDelSide=2092
return
!
20 iold=dRef(i1)
iend=dNeib(iold)%next
i=iend
!
30   if(dNeib(i)%neib.EQ.i0) then
 call DelNeib(i1,iold,i)
 goto 40
end if
iold=i
i=dNeib(iold)%next
if(i.ne.iend) goto 30
!
iDelSide=2093
return
!
40 iold=dRef(i2)
neib1=dNeib(iold)%neib
iend=dNeib(iold)%next
i=iend
!
50   neib0=neib1
neib1=dNeib(i)%neib
if(neib0.EQ.i0 .and. neib1.EQ.i1) then
 ii=dPusto
 dPusto=dNeib(ii)%next
 if(dPusto.EQ.0) then
  iDelSide=2093
  return
 end if
 nFikt=nFikt+1
 dNeib(iold)%next=ii
 dNeib(ii)%next=i
 dNeib(ii)%neib=nFikt
 rq=2*((x(i2)-x0)*dy-(y(i2)-y0)*dx)*d2
 x(nFikt)=x(i2)-rq*dy
 y(nFikt)=y(i2)+rq*dx
 dType(nFikt-n)=i0
 goto 60
end if
iold=i
i=dNeib(iold)%next
if(i.ne.iend) goto 50
!
iDelSide=2094
return
!
60 iDelSide=0
end
!
!**********************************************************************
!
integer function iTestTri(it1,i1,ig1)
use Mesh_data
integer it1,i1,ig1
double precision x0,y0,x1,y1,dx,dy,d2
common/dDir2/ x0,y0,x1,y1,dx,dy,d2
double precision x2,y2,STri
integer iv,it
integer i,rc
integer TriList(3,500),nList,iList
!!integer TriList(3,100),nList,iList
integer it0,i0,ig
!
integer iDelSide,iGetVerN! Functions
!
nList=1
TriList(1,1)=it1
TriList(2,1)=i1
TriList(3,1)=ig1
iList=1
!
10 it0=TriList(1,iList)
i0 =TriList(2,iList)
ig =TriList(3,iList)
!
iv=tri(it0)%v(i0)
x2=xv(it0)
y2=yv(it0)
STri=dx*(y0-y2)-dy*(x0-x2)
if(STri.LT.0.) then
 !
 rc=iDelSide(it0,i0)
 if(rc.ne.0) then
  iTestTri=rc
  return
 end if
 !
 it=tri(it0)%t(LEFT(i0))
 i=iGetVerN(it,ig)
 nList=nList+1
 TriList(1,nList)=it
 TriList(2,nList)=LEFT(i)
 TriList(3,nList)=ig
 !
 it=tri(it0)%t(RIGHT(i0))
 i=iGetVerN(it,iv)
 nList=nList+1
 TriList(1,nList)=it
 TriList(2,nList)=LEFT(i)
 TriList(3,nList)=iv
 !
end if
!
if(iList.lt.nList) then
 iList=iList+1
 goto 10
end if
!
iTestTri=0
end
!
!**********************************************************************
!
function iSgCorr()
use Mesh_data
double precision x0,y0,x1,y1,dx,dy,d2
common/dDir2/ x0,y0,x1,y1,dx,dy,d2
integer i,j
integer ig,ig1,it,rc
!
integer iGetVerN,iTestTri,iSgCorr   ! Functions
!
nFikt=n+3
!
do ig=1,ng
 ig1=ig+1
 if(ig.EQ.ng) ig1=1
 x0=x(ig)
 x1=x(ig1)
 y0=y(ig)
 y1=y(ig1)
 dx=x0-x1
 dy=y0-y1
 d2=1./(dx**2+dy**2)
 !
 it=ftri(ig)
 !
 10     i=iGetVerN(it,ig)
 j=RIGHT(i)
 !
 if(tri(it)%v(j).EQ.ig1) goto 20
 it=tri(it)%t(j)
 !
 goto 10
 !
 20   rc=iTestTri(it,int(LEFT(i)),ig)
 if(rc.ne.0) then
  iSgCorr=rc
  return
 end if
end do
!
iSgCorr=0
end
