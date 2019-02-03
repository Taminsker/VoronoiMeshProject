module our_module

  use mesh_data

  implicit none

contains
function areaTri(p1, p2, p3)
  real*8, dimension(0:1),  intent(inout) :: p1,p2,p3
  real*8 :: areaTri
  areaTri = (0.5) * abs((p2(0)-p1(0))*(p3(1)-p1(1))- (p3(0)-p1(0))*(p2(1)-p1(1)))
  print*,' areaTri = ', areaTri

end function areaTri

subroutine isoBArea(isoB, cnum,Mesh,XYp)
  implicit none
  !-------------
  ! Parameters
  type(Mesh_struct),     intent(in) :: Mesh
  real*8, dimension(:,:),  intent(in) :: XYp
  real*8, dimension(0:1)            :: sumForIso, p1, p2
  real*8                            :: sumForVol
  character(10), intent(inout)      :: isoB
  character(10)                    :: NextFile

  character(30) ,intent(inout)      :: cnum
  integer :: i,j
  real*8 :: xx,yy
  !------------
  NextFile = 'isoB.'
  isoB = trim(NextFile)//trim(cnum)
  isoB = trim(isoB)

  open(25, file=isoB)

  xx = 0.9d0
  xx = 1.0d0
  yy = 1.0d0 - xx
  do j = 1,Mesh%nc
    sumForIso = 0
    sumForVol = 0

     do i = 1,Mesh%c_l(j)
       sumForIso(0) = sumForIso(0) + xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1)
       sumForIso(1) = sumForIso(1) + xx* yn( Mesh%cell_list(j,i) ) + yy * XYp(j,2)
     end do

     sumForIso(0) = sumForIso(0) / Mesh%c_l(j)
     sumForIso(1) = sumForIso(1) / Mesh%c_l(j)

     do i = 1,Mesh%c_l(j)
       if (j.NE.Mesh%nc) then
       p1(0)= xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1)
       p1(1)= xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,2)

       p2(0)= xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j+1,1)
       p2(1)= xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j+1,2)
       sumForVol = sumForVol + areaTri(sumForIso, p1, p2)
     end if
     end do

     p1(0)= xx* xn( Mesh%cell_list(j,1) ) + yy * XYp(j,1)
     p1(1)= xx* xn( Mesh%cell_list(j,1) ) + yy * XYp(j,2)

     p2(0)= xx* xn( Mesh%cell_list(j,Mesh%c_l(j)) ) + yy * XYp(j,1)
     p2(1)= xx* xn( Mesh%cell_list(j,Mesh%c_l(j)) ) + yy * XYp(j,2)
     sumForVol = sumForVol + areaTri(sumForIso, p1, p2)

     write(25,*) sumForIso(0), sumForIso(1), sumForVol
  end do
  close(25)
end subroutine isoBArea

end module our_module
