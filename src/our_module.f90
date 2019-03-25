module our_module

  use mesh_data

  implicit none

contains
  function areaTri(p1, p2, p3)
    real*8,  dimension(0:1), intent(inout) :: p1, p2, p3
    real :: areaTri
    areaTri = (0.5d0) * abs((p2(0)-p1(0))*(p3(1)-p1(1))- (p3(0)-p1(0))*(p2(1)-p1(1)))
  end function areaTri

  subroutine isoBArea(isoB, cnum, Mesh, XYp)
    implicit none
    !-------------
    ! Parameters
    type(Mesh_struct),     intent(in) :: Mesh
    real*8, dimension(:,:),  intent(in) :: XYp
    real*8, dimension(0:1)            :: sumForIso, p1, p2, p3
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
      end do

      sumForIso(0) = sumForIso(0) / Mesh%c_l(j)
      sumForIso(1) = sumForIso(1) / Mesh%c_l(j)
      ! write(25,*) sumForIso
      ! write(25,*), ' '

      do i = 1,Mesh%c_l(j)
        if (i.NE.Mesh%c_l(j)) then
          p1 = sumForIso
          p2(0) = xx* xn( Mesh%cell_list(j,i) ) + yy * XYp(j,1)
          p2(1) = xx* yn( Mesh%cell_list(j,i) ) + yy * XYp(j,2)

          p3(0) = xx* xn( Mesh%cell_list(j,i+1) ) + yy * XYp(j,1)
          p3(1) = xx* yn( Mesh%cell_list(j,i+1) ) + yy * XYp(j,2)

          sumForVol = sumForVol + areaTri(p1, p2, p3)
        end if
      end do
      ! write(25,*), ' '

      p1 = sumForIso
      p2(0) = xx* xn( Mesh%cell_list(j,1) ) + yy * XYp(j,1)
      p2(1) = xx* yn( Mesh%cell_list(j,1) ) + yy * XYp(j,2)
      p3(0) = xx* xn( Mesh%cell_list(j,Mesh%c_l(j)) ) + yy * XYp(j,1)
      p3(1) = xx* yn( Mesh%cell_list(j,Mesh%c_l(j)) ) + yy * XYp(j,2)

      sumForVol = sumForVol + areaTri(p1, p2, p3)
      volumeTotal = volumeTotal + sumForVol;
      write(25,*) sumForIso!, sumForVol
    end do
    print*, 'volume total = ', volumeTotal
    close(25)
  end subroutine isoBArea


  function potential(t,t0,t1,a,b)
    implicit none
    real*8, intent(in) :: a, b
    integer, intent(in) :: t0,t1,t
    real*8 :: potential

    if (t < t0) then
      potential = a
    elseif (t < t1) then
      potential = ((a - b) / (t0 - t1) * t) + a - (a - b) / (t0 - t1) * t0
    else
      potential = b
    end if
  end function potential

  subroutine zone(Mesh, XYp, x0, y0, d0, x1, y1, d1)
    implicit none

    type(Mesh_struct),     intent(inout) :: Mesh
    real*8, intent(in) :: x0, y0, d0, d1, x1, y1
    real*8,  dimension(:,:) :: XYp

    integer :: i

    do i = 5, 105
      call random_number(XYp(i,1))
      call random_number(XYp(i,2))
      XYp(i,1) = XYp(i,1) * d0 + x0
      XYp(i,2) = XYp(i,2) * d0 + y0
    end do

    do i = 106, 206
      call random_number(XYp(i,1))
      call random_number(XYp(i,2))
      XYp(i,1) = XYp(i,1) * d1 + x1
      XYp(i,2) = XYp(i,2) * d1 + y1
    end do
  end subroutine

  ! subroutine temporaire(Mesh)
  !   implicit none
  !   type(Mesh_struct),     intent(inout) :: Mesh
  !   open(26,file="neighbor.txt")
  !   write(26,*) Mesh%neighbor
  !   close(26)
  ! end subroutine


  subroutine odour_transmission(Mesh,time, W_c, W_p)
    implicit none
    type(Mesh_struct), intent(in) :: Mesh
    real*8,  dimension(:), intent(inout) :: W_c,W_p
    integer, intent (in) :: time
    ! real*8 :: k
    integer :: i,j, icycle2_the_return_MOMY, k

    ! initialisation step

    W_c(1:4) = 0
    W_c(5:104) = potential(time, 133, 233, 0.25d0, 0.75d0)
    W_c(105:204) = potential(time, 100, 350, 3d0, 0d0)
    W_c(205:Mesh%nc)=0

    ! beginning of the odour transmition

    do icycle2_the_return_MOMY = 1,Mesh%nc
      do i = 1,Mesh%nn
        W_p(i) = 0.0d0
        do j = 1,Mesh%n_l(i)  ! boucle sur les noeuds
          k = Mesh%node_list(i,j)  ! numero de la j-ieme cellule autours du noeuds
          W_p(i) = W_p(i) + W_c(k)
        end do
        W_p(i) = W_p(i)/dble(Mesh%n_l(i)) ! dble convert an int into a double
      end do

      do j = 1,Mesh%nc
        W_c(j) = 0
        do i= 1,Mesh%c_l(j)
          k = Mesh%cell_list(j,i)
          W_c(j) = W_c(j) + W_p(k)
        end do
        W_c(j) = W_c(j)/Mesh%c_l(j)
      end do
    end do

    W_c(1:4) = 0
    W_c(5:104) = potential(time, 133, 233, 0.25d0, 0.75d0)
    W_c(105:204) = potential(time, 100, 350, 1d0, 0d0)

  end subroutine
end module our_module
