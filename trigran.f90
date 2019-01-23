      subroutine Bound
use mesh_data

      integer k1,k2,i

      k2=ng

      !print*,'BOUND',ng

      do 100 i=1,ng
        k1=k2
        k2=i
        call TriGran(k1,k2)
  100 continue

      do i=1,n
        ftri(i)=0
      end do
      call GetFirst

      end

!**************************************************************

integer function iGetVerN(it,iv)
use mesh_data
      integer it,iv

      if(tri(it)%v(0).eq.iv) then
        iGetVerN=0
        return
      endif

      if(tri(it)%v(1).eq.iv) then
        iGetVerN=1
        return
      endif

      iGetVerN=2

      return
      end

!**************************************************************

      subroutine TriGran(p1,p2)
use mesh_data
      integer p1,p2
      integer it,it0,j,jr,vr,vl,vm
      integer kt
      double precision yr,yl
      double precision dist,dx,dy
      logical first
      integer iGetVerN             ! Function

      integer triList(0:500)
      integer vmList(0:500)
      integer vpList(0:500)
      integer tmList(0:500)
      integer tpList(0:500)
      integer nTriList,nVmList,nVpList
      double precision cosFi,sinFi,x_0,y_0
      integer p2_0,it_0

      common/dGran/ triList,vmList,vpList,tmList,tpList, &
                    nTriList,nVmList,nVpList,            &
                    cosFi,sinFi,x_0,y_0,p2_0,it_0

      nTriList=0
      nVmList=1
      nVpList=1
      vmList(0)=p1
      vpList(0)=p1

      x_0=x(p1)
      y_0=y(p1)
      dx=x(p2)-x(p1)
      dy=y(p2)-y(p1)
      dist=dsqrt(dx*dx+dy*dy)
      cosFi=dx/dist
      sinFi=dy/dist
      yr=0
      yl=0

!
! 1. Find starting triangle with vertex p1
!

      it=ftri(p1)
      j=iGetVerN(it,p1)
      if(tri(it)%v(j).ne.p1) then
        if(p2_0.eq.p1) then
          it=it_0
        else
          do it=1,nt-1
            if(tri(it)%v(0).eq.p1 .or.     &
               tri(it)%v(1).eq.p1 .or.     &
               tri(it)%v(2).eq.p1) goto 100
          end do
        end if
      end if
  100 continue

      p2_0=p2

!
! 2. Try to find triangle with vertixe p1 in wich lies segment (p1,p2)
!

      first=.true.

  200 continue

        j=iGetVerN(it,p1)
        jr=RIGHT(j)
        vr=tri(it)%v(jr)
        vl=tri(it)%v(LEFT(j))
        if(vr.eq.p2 .or. vl.eq.p2) then
          it_0=it
          return
        end if
        if(first) then
          yr=-(x(vr)-x_0)*sinFi+(y(vr)-y_0)*cosFi
        else
          yr=yl
        end if
        yl=-(x(vl)-x_0)*sinFi+(y(vl)-y_0)*cosFi
        it0=it

        it=tri(it)%t(jr)

      if(yr.ge.0. .or. yl.lt.0.) goto 200

      if(yl.eq.0.) return                       ! Point on the line

!
! 3. Go through triangles along segment (p1,p2)
!

      triList(nTriList)=it0
      nTriList=nTriList+1
      it=tri(it0)%t(j)
      tmList(nVmList)=tri(it0)%t(LEFT(j))
      vmList(nVmList)=vr
      vm=vr
      tpList(nVpList-1)=tri(it0)%t(jr)
      vpList(nVpList)=vl
      nVmList=nVmList+1
      nVpList=nVpList+1

  300 continue
        triList(nTriList)=it
        nTriList=nTriList+1
        j=iGetVerN(it,vm)
        vr=tri(it)%v(RIGHT(j))
        if(vr.eq.p2) goto 400
        yr=-(x(vr)-x_0)*sinFi+(y(vr)-y_0)*cosFi
        if(yr.gt.0) then
          tpList(nVpList-1)=tri(it)%t(j)
          vpList(nVpList)=vr
          nVpList=nVpList+1
          it=tri(it)%t(LEFT(j))
        else
          tmList(nVmList)=tri(it)%t(LEFT(j))
          vmList(nVmList)=vr
          vm=vr
          nVmList=nVmList+1
          it=tri(it)%t(j)
        end if
      goto 300

  400 continue

      tmList(nVmList)=tri(it)%t(LEFT(j))
      tpList(nVpList-1)=tri(it)%t(j)
      vmList(nVmList)=p2
      nVmList=nVmList+1
      vpList(nVpList)=p2
      nVpList=nVpList+1

!
!  Change triangulations under and over segment (p1,p2)
!

      call ReTriang(p1,vmList,tmList,nVmList,1)
      call ReTriang(p2,vpList,tpList,nVpList,0)

!
!  Connect two last triangles
!

      it=tpList(0)
      kt=tmList(nVmList-1)
      tri(it)%t(LEFT(iGetVerN(it,p1)))=kt
      tri(kt)%t(LEFT(iGetVerN(kt,p2)))=it

      it_0=it

      return
      end

!**************************************************************

      subroutine ReTriang(p1,v_List,t_List,n_List,plus)
use mesh_data
      integer p1,v_List(0:500),t_list(0:500),n_List,plus
      logical minus
      integer it,j,vr,vl,vm,v
      integer tt,tr,tl,k,kr,kl
      integer i,ii,kt
      double precision dist
      integer iGetVerN

      integer triList(0:500)
      integer vmList(0:500)
      integer vpList(0:500)
      integer tmList(0:500)
      integer tpList(0:500)
      integer nTriList,nVmList,nVpList
      double precision cosFi,sinFi,x_0,y_0
      double precision Stri,Dist2
      integer p2_0,it_0

      common/dGran/ triList,vmList,vpList,tmList,tpList,  &
                    nTriList,nVmList,nVpList,             &
                    cosFi,sinFi,x_0,y_0,p2_0,it_0

      kt=n_List-2
      minus=(plus.eq.0)

  100 if(kt.le.0) return

        vr=p1
        v=-1
        kr=0
        do 300 ii=1,n_List-1
          i=ii
          if(minus) i=n_List-ii-1
          vm=v_List(i)
          if(vm.lt.0) goto 300
          vl=v
          v=vr
          vr=vm
          tl=tt
          tt=tr
          tr=t_List(i)
          kl=k
          k=kr
          kr=i
          if(vl.gt.0 .and. Stri(vl,v,vr).gt.0) then
            it=triList(nTriList-1)
            call GetV(vl,v,vr,it)
            dist=Dist2(it,v)

            do 200 j=0,n_List-1
              vm=v_List(j)
              if(vm.lt.0) goto 200
              if(vl.eq.vm .or. v.eq.vm .or. vr.eq.vm) goto 200
              if(Dist2(it,vm).lt.dist) goto 300
  200       continue

            nTriList=nTriList-1
            kt=kt-1
            tri(it)%v(0)=vl
            tri(it)%v(1)=v
            tri(it)%v(2)=vr
            tri(it)%t(0)=tr
            tri(it)%t(1)=-1
            tri(it)%t(2)=tt

            ftri(vl)=it
            ftri(v )=it
            ftri(vr)=it

            j=iGetVerN(tt,v)
            tri(tt)%t(LEFT(j))=it
            j=iGetVerN(tr,v)
            tri(tr)%t(RIGHT(j))=it
            v=-v
            v_List(k)=v
            tr=it
            t_List(kr)=it

          end if

  300   continue

      goto 100
      end

!**************************************************************

      double precision function Dist2(tt,v)
use mesh_data
      integer tt,v

      Dist2=(xv(tt)-x(v))**2+(yv(tt)-y(v))**2

      return
      end

!**************************************************************

      double precision function Stri(i1,i2,i3)
use mesh_data
      integer i1,i2,i3
      double precision x1,x2,x3,y1,y2,y3

      x1=x(i1)
      y1=y(i1)
      x2=x(i2)
      y2=y(i2)
      x3=x(i3)
      y3=y(i3)
      Stri=(x1-x2)*(y1-y3)-(y1-y2)*(x1-x3)
      return
      end

