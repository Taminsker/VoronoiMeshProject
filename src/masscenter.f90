module masscenter

  use mesh_data
  use our_module

  implicit none

contains

function triMassCenter(p1, p2, p3)
  real*8, dimension(0:1), intent(in) :: p1, p2, p3
  real*8, dimension(0:1) :: triMassCenter

! p1 --> p2
  triMassCenter(0) = 1d0/6d0 * (p1(0) * p1(0) + p1(0) * p2(0) + p2(0) * p2(0)) * (p2(1) - p1(1))
  triMassCenter(1) = -1d0/6d0 * (p1(1) * p1(1) + p1(1) * p2(1) + p2(1) * p2(1)) * (p2(0) - p1(0))
  print*, 'triMassCenter 1 = ', triMassCenter
! p2 --> p3
  triMassCenter(0) = triMassCenter(0) + 1d0/6d0 * (p2(0) * p2(0) + p2(0) * p3(0) + p3(0) * p3(0)) * (p3(1) - p2(1))
  triMassCenter(1) = triMassCenter(1) - 1d0/6d0 * (p2(1) * p2(1) + p2(1) * p3(1) + p3(1) * p3(1)) * (p3(0) - p2(0))
  print*, 'triMassCenter 2 = ', triMassCenter

! p3 --> p1
  triMassCenter(0) = triMassCenter(0) + 1d0/6d0 * (p3(0) * p3(0) + p3(0) * p1(0) + p1(0) * p1(0)) * (p1(1) - p3(1))
  triMassCenter(1) = triMassCenter(1) - 1d0/6d0 * (p3(1) * p3(1) + p3(1) * p1(1) + p1(1) * p1(1)) * (p1(0) - p3(0))
  print*, 'triMassCenter 3 = ', triMassCenter


end function triMassCenter

subroutine cellMassCenter(massCenter, cnum, Mesh, XYp)
  implicit none
  !-------------
  ! Parameters
  type(Mesh_struct),     intent(in) :: Mesh
  real*8, dimension(:,:),intent(in) :: XYp
  real*8, dimension(0:1)            :: sumForMass, p1, p2, p3
  real*8                            :: sumForVol
  real*8                            :: volumeTotal
  character(10), intent(inout)      :: massCenter
  character(10)                     :: NextFile

  character(30) ,intent(inout)      :: cnum
  integer :: i,j
  real*8 :: xx,yy
  !------------
  volumeTotal = 0d0
  NextFile = 'massCenter.'
  massCenter = trim(NextFile)//trim(cnum)
  massCenter = trim(massCenter)

  open(25, file=massCenter)
  open(26, file='triMassCenter')

  ! do j = 1,Mesh%nc
  do j = 5,5

    sumForMass = 0d0
    sumForVol = 0d0

    ! do i = 1,Mesh%c_l(j)
    !   sumForMass(0) = sumForMass(0) + xn( Mesh%cell_list(j,i) )
    !   sumForMass(1) = sumForMass(1) + yn( Mesh%cell_list(j,i) )
    ! end do
    !
    ! sumForMass(0) = sumForMass(0) / Mesh%c_l(j)
    ! sumForMass(1) = sumForMass(1) / Mesh%c_l(j)
    ! ! write(25,*) sumForMass
    ! ! write(25,*), ' '

    do i = 1,Mesh%c_l(j)
      if (i.NE.Mesh%c_l(j)) then
        p1(0) = XYp(j,1)
        p1(1) = XYp(j,2)

        p2(0) = xn( Mesh%cell_list(j,i) )
        p2(1) = yn( Mesh%cell_list(j,i) )

        p3(0) = xn( Mesh%cell_list(j,i+1) )
        p3(1) = yn( Mesh%cell_list(j,i+1) )

        sumForMass = sumForMass + areaTri(p1, p2, p3) * triMassCenter(p1, p2, p3)
        ! write(26,*) triMassCenter(p1, p2, p3)
        print*, 'p1 = ', p1
        print*, 'p2 = ', p2
        print*, 'p3 = ', p3, ' '
        write(26,*) p1
        write(26,*) p2
        write(26,*) p3


        sumForVol = sumForVol + areaTri(p1, p2, p3)
      end if
      ! print*,i,'COUCOU1'
    end do
    ! write(25,*), ' '

    p1(0) = XYp(j,1)
    p1(1) = XYp(j,2)
    ! print*,j,'COUCOU - 3'

    p3(0) = xn( Mesh%cell_list(j,1) )
    p3(1) = yn( Mesh%cell_list(j,1) )
    ! print*,'COUCOU - 5'

    p2(0) = xn( Mesh%cell_list(j,Mesh%c_l(j)) )
    p2(1) = yn( Mesh%cell_list(j,Mesh%c_l(j)) )

    sumForMass = sumForMass + areaTri(p1, p2, p3) * triMassCenter(p1, p2, p3)
    sumForVol = sumForVol + areaTri(p1, p2, p3)

    volumeTotal = volumeTotal + sumForVol
    sumForMass = sumForMass / volumeTotal
    print*, 'mass center j =', j, ' : ', sumForMass, 'areaTri = ', areaTri(p1, p2, p3)
    write(25,*) sumForMass!, sumForVol
  end do
  print*, 'volume total = ', volumeTotal
  close(25)
  close(26)
end subroutine cellMassCenter
end module masscenter
