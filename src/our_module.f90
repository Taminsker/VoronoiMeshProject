module our_module

  use mesh_data

  implicit none

contains
  function areaTri(p1_0, p1_1, p2_0, p2_1, p3_0, p3_1)
    real*8,  intent(inout) :: p1_0, p1_1, p2_0, p2_1, p3_0, p3_1
    real :: areaTri
    areaTri = (0.5d0) * abs((p2_0-p1_0)*(p3_1-p1_1)- (p3_0-p1_0)*(p2_1-p1_1))
    ! print*, 'p1 =', p1_0, p1_1
    ! print*, 'p2 =', p2_0, p2_1
    ! print*, 'p3 =', p3_0, p3_1
    ! print*,' areaTri = ', areaTri
    ! print*, ' '

  end function areaTri

  subroutine isoBArea(isoB, cnum,Mesh,XYp)
    implicit none
    !-------------
    ! Parameters
    type(Mesh_struct),     intent(in) :: Mesh
    real*8, dimension(:,:),  intent(in) :: XYp
    real*8, dimension(0:1)            :: sumForIso
    real*8                            :: p1_0, p1_1, p2_0, p2_1, p3_0, p3_1
    real*8                            :: sumForVol
    real*8                            :: volumeTotal
    character(10), intent(inout)      :: isoB
    character(10)                    :: NextFile

    character(30) ,intent(inout)      :: cnum
    integer :: i,j
    real*8 :: xx,yy
    !------------
    volumeTotal = 0d0
    NextFile = 'isoB.'
    isoB = trim(NextFile)//trim(cnum)
    isoB = trim(isoB)

    open(25, file=isoB)

    xx = 0.9d0
    xx = 1.0d0
    yy = 1.0d0 - xx
    do j = 1,Mesh%nc
      sumForIso = 0d0
      sumForVol = 0d0

      do i = 1,Mesh%c_l(j)
        sumForIso(0) = sumForIso(0) + xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1)
        sumForIso(1) = sumForIso(1) + xx* yn( Mesh%cell_list(j,i) ) + yy * XYp(j,2)
        ! write(25,*), xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1), xx* yn( Mesh%cell_list(j,i) ) + yy * XYp(j,2)
      end do

      sumForIso(0) = sumForIso(0) / Mesh%c_l(j)
      sumForIso(1) = sumForIso(1) / Mesh%c_l(j)
      write(25,*), sumForIso(0), sumForIso(1)
      ! write(25,*), ' '

      do i = 1,Mesh%c_l(j)
        if (i.NE.Mesh%c_l(j)) then
          p1_0 = sumForIso(0)
          p1_1 = sumForIso(1)

          p2_0 = xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1)
          p2_1 = xx* yn( Mesh%cell_list(j,i) ) + yy * XYp(j,2)

          p3_0 = xx* xn( Mesh%cell_list(j,i+1) ) + yy * XYp(j,1)
          p3_1 = xx* yn( Mesh%cell_list(j,i+1) ) + yy * XYp(j,2)

          sumForVol = sumForVol + areaTri(p1_0, p1_1, p2_0, p2_1, p3_0, p3_1)
        end if
      end do
      ! write(25,*), ' '

      p1_0 = sumForIso(0)
      p1_1 = sumForIso(1)
      p2_0 = xx* xn( Mesh%cell_list(j,1) ) + yy * XYp(j,1)
      p2_1 = xx* yn( Mesh%cell_list(j,1) ) + yy * XYp(j,2)
      p3_0 = xx* xn( Mesh%cell_list(j,Mesh%c_l(j)) ) + yy * XYp(j,1)
      p3_1 = xx* yn( Mesh%cell_list(j,Mesh%c_l(j)) ) + yy * XYp(j,2)

      sumForVol = sumForVol + areaTri(p1_0, p1_1, p2_0, p2_1, p3_0, p3_1)
      volumeTotal = volumeTotal + sumForVol;
      ! write(25,*) sumForIso(0), sumForIso(1)!, sumForVol
    end do
    print*, 'volume total = ', volumeTotal
    close(25)
  end subroutine isoBArea

end module our_module
