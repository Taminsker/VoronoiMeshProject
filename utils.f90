! Get the number of the neighbours for the point i

      integer function iGetNeibN(i)
use mesh_data 
!
      integer i
      integer nn0,ind,ibeg

      ind=dRef(i)
      ibeg=ind
      nn0=0

  100 continue

        ind=dNeib(ind)%next
        nn0=nn0+1
        if(ind.ne.ibeg) goto 100

      iGetNeibN=nn0
      end

!------------------------------------------------------------------------

! Store neighbours of the point i in the array
! Additionaly store first neighbour in the last position:
!
!         N1 N2 N3 ... Nn N1

      subroutine GetNbAr(i,array,na)
use mesh_data 
!
      integer i,na,array(na)
!
! i     - Number of point (in)
! array - Array of numbers of the neighbours (out)
! na    - Length of the array (in)
!
      integer ip,ind,ibeg

      ind=dRef(i)
      ibeg=ind
      ip=1
      if(ip.gt.na) return
      array(ip)=dNeib(ind)%neib
      !print*,ip,ind,array(ip)

  100 continue

        ind=dNeib(ind)%next
        ip=ip+1
        if(ip.gt.na) return
        array(ip)=dNeib(ind)%neib

        !print*,ip,ind,array(ip)

        if(ind.ne.ibeg) goto 100

      end

!------------------------------------------------------------------------

! Store vertices of the cell I in the arrays X and Y:
! Additionaly store first vertice in the last position:
!
!         X1 X2 X3 ... Xn X1
!         Y1 Y2 Y3 ... Yn Y1
!
! (X1,Y1) is a vertice for the triangle (I,N1,N2) - see subroutine
! GetNbAr

      subroutine GetVAr(i,coordX,coordY,na)
use mesh_data 
!
      integer i,na
      double precision coordX(na),coordY(na)
!
! i      - Number of point (in)
! coordX - Array of coordinates X of the cell's vertices (out)
! coordY - Array of coordinates Y of the cell's vertices (out)
! na     - Length of the arrays (in)
!
      integer ip,ind,ibeg

      ind=dRef(i)
      ibeg=ind
      ip=1
      if(ip.gt.na) return
      coordX(ip)=xVert(ind)
      coordY(ip)=yVert(ind)


  100 continue

        ip=ip+1
        ind=dNeib(ind)%next
        if(ip.gt.na) return
        coordX(ip)=xVert(ind)
        coordY(ip)=yVert(ind)
        if(ind.ne.ibeg) goto 100

      end

