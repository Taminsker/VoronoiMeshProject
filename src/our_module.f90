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
    real*8 :: a, W_t
    integer :: i,j, icycle2_the_return_MOMY, k
    integer, dimension(1:Mesh%nc) :: cent
    integer, dimension(1:Mesh%nn) :: nodes
    logical :: anotherRound
    anotherRound = .TRUE.


    ! do i = 1,Mesh%nc
    !   a = (Mesh%X_c(i)-0.3d0)**2 + (Mesh%Y_c(i)-0.3d0)**2
    !   W_c(i) = potential(time, 133, 160, 0.25d0, 0.75d0) * exp(-a/0.5d0)
    !
    !   a = (Mesh%X_c(i)-0.8d0)**2 + (Mesh%Y_c(i)-0.8d0)**2
    !   W_c(i) = W_c(i) + potential(time, 133, 160, 0.75d0, 0.25d0) * exp(-a/0.5d0)
    !   ! if (W_c(i) <= 1e-10) then
    !   !   print*, i
    !   ! end if
    ! end do
    !
    ! do i = 1,Mesh%nn
    !   a = (Mesh%X_n_n(i)-0.3d0)**2 + (Mesh%Y_n_n(i)-0.3d0)**2
    !   W_p(i) = potential(time, 133, 160, 0.25d0, 0.75d0) * exp(-a/0.5d0)
    !
    !   a = (Mesh%X_n_n(i)-0.8d0)**2 + (Mesh%Y_n_n(i)-0.8d0)**2
    !   W_p(i) = W_p(i) + potential(time, 133, 160, 0.75d0, 0.25d0) * exp(-a/0.5d0)
    !   ! if (W_p(i) <= 1e-10) then
    !   !   print*, i
    !   ! end if
    ! end do
    ! initialisation step
    W_c = 0d0
    W_p = 0d0



    ! beginning of the odour transmition

    cent = 0
    nodes = 0

    do while (anotherRound .EQV. .TRUE.)
    ! do icycle2_the_return_MOMY = 1,nint(Mesh%nc/2d0)+1   ! Modifier le nb d'itérations pour passer les grands nombres de générateurs

      if (time >= t0_source1) then
        W_c(5:104) = potential(time, t0_source1, t1_source1, a_source1, b_source1)
      end if
      if (time <= t1_source2) then
        W_c(105:204) =potential(time, t0_source2, t1_source2, a_source2, b_source2)
      end if

      do i = 1,Mesh%nn
        W_t = 0.0d0
        do j = 1,Mesh%n_l(i)  ! boucle sur les noeuds
          k = Mesh%node_list(i,j)  ! numero de la j-ieme cellule autours du noeuds
          W_t = W_t + W_c(k)
        end do

        if (.NOT. W_t == 0d0) then
          nodes(i) = 1
        end if

        ! W_p(i) = (W_p(i) + W_t)/dble(Mesh%n_l(i)) ! dble convert an int into a double
        W_p(i) = W_t/dble(Mesh%n_l(i)) ! dble convert an int into a double

      end do

      do j = 1,Mesh%nc
        W_t = 0d0
        do i= 1,Mesh%c_l(j)
          k = Mesh%cell_list(j,i)
          W_t = W_t + W_p(k)
        end do

        if (.NOT. W_t == 0d0) then
          cent(j) = 1
        end if

        ! W_c(j) = (W_c(j) + W_t)/dble(Mesh%c_l(j))
        W_c(j) = W_t/dble(Mesh%c_l(j))

      end do

      anotherRound = .FALSE.
      do i = 1,Mesh%nn
        if (nodes(i) == 0) then
          anotherRound = .TRUE.
          EXIT
        end if
      end do
      do i = 1,Mesh%nc
        if (cent(i) == 0) then
          anotherRound = .TRUE.
          EXIT
        end if
      end do

    end do

    ! W_c(1:4) = 0d0
    ! if (time >= 100) then
    !   W_c(5:104) = potential(time, 100, 160, 0.25d0, 0.9d0)
    ! end if
    ! if (time <= 210) then
    !   W_c(105:204) = potential(time, 133, 210, 0.9d0, 0.25d0)
    ! end if
    !
    ! print*, 'W_c(205:Mesh%nc) = ', W_c(205:Mesh%nc)
    ! ! print*, 'W_p() = ', W_p
  end subroutine



  subroutine elephantOdourContinuous(Mesh, W_c, W_p, time)
    implicit none
    type(Mesh_struct), intent(inout) :: Mesh
    real*8,  dimension(:), intent(in) :: W_c,W_p
    real*8 :: x, y, xprime, yprime
    integer :: i, j
    real*8, dimension(2) :: newGen, Cn, temp
    integer, intent (in) :: time
    ! real*8 :: maxOdourOnCell

    ! temp = (/ potential(time, 0, 466, 0.25d0, 0.75d0), potential(time, 200, 600, 1d0, 0d0)/)


    do i = 205, Mesh%nc
      newGen = 0d0

      do j = 1,Mesh%c_l(i)
        x = sqrt((xn(Mesh%cell_list(i,j)) - x_lac)**2 + (yn(Mesh%cell_list(i,j))  - y_lac)**2)
        if (x > rayon) then
          Cn = (/xn(Mesh%cell_list(i,j)) - Mesh%X_c(i), yn(Mesh%cell_list(i,j)) - Mesh%Y_c(i)/)
          Cn = Cn / (sqrt(Cn(1)**2 + Cn(2)**2))
          newGen = newGen + Cn * DMAX1(0d0, W_p(Mesh%cell_list(i,j)) - W_c(i))
          ! newGen = newGen + Cn * (W_p(Mesh%cell_list(i,j)) - W_c(i))
        end if
      end do

      if (.NOT. sqrt(newGen(1)**2 + newGen(2)**2) == 0d0) then
      newGen = newGen / (100 / xmax * sqrt(newGen(1)**2 + newGen(2)**2)) + (/ Mesh%X_c(i), Mesh%Y_c(i) /)
      else
        newGen = (/ Mesh%X_c(i), Mesh%Y_c(i) /)
      end if

      x = sqrt((newGen(1) - x_lac)**2 + (newGen(2) - y_lac)**2)

      if (x <= rayon) then
        call random_number(y)
        y = nint(y * 10)
        newGen(1) = x_lac + (rayon) * (newGen(1) - x_lac) / x + y * 1e-3
        newGen(2) = y_lac + (rayon) * (newGen(2) - y_lac) / x + y * 1e-3
      end if

      Mesh%X_c(i) = newGen(1)
      Mesh%Y_c(i) = newGen(2)

      ! if ((.NOT. Mesh%X_c(i) == DMAX1(0d0 + 1e-8, DMIN1(1d0 - 1e-8 , newGen(1)))) &
      ! .OR. (.NOT. Mesh%Y_c(i) == DMAX1(0d0 + 1e-8, DMIN1(1d0 - 1e-8 , newGen(2))))) then
      !   print*, Mesh%X_c(i), Mesh%Y_c(i)
      ! end if
      if (isnan(newGen(1))) stop '"newGen(1)" is a NaN'
      if (isnan(newGen(2))) stop '"newGen(2)" is a NaN'

      if (depassement == 0) then
        Mesh%X_c(i) = DMAX1(xmin + 1e-8, DMIN1(xmax - 1e-8 , newGen(1)))
        Mesh%Y_c(i) = DMAX1(ymin + 1e-8, DMIN1(ymax - 1e-8, newGen(2)))
      end if




      ! print*, 'after :'
      ! print*, "x = ", Mesh%X_c(i)
      ! print*, "y = ", Mesh%Y_c(i)
    end do

  end subroutine

  subroutine creat_script3(icycle, Mesh, XYp)

    implicit none
    type(Mesh_struct), intent(inout) :: Mesh
    integer, intent(in) :: icycle
    integer :: p, i, j
    real*8,  dimension(:,:),   intent(in) :: XYp
    real*8 :: x, y, xprime, yprime
    character(30)       :: cnum, filename

    write(cnum,*) icycle + 10000
    cnum     = adjustl(cnum)

    ! ----------------------------------------
    filename  = 'trees.'
    filename = trim(filename)//trim(cnum)

    open(21,file=filename)
    do p = 1,204
       write(21,*) XYp(p,1),XYp(p,2)
    end do
    close(21)

    ! ----------------------------------------
    filename  = 'elephants.'
    filename = trim(filename)//trim(cnum)

    open(21,file=filename)
    do p = 205,Mesh%nc
       write(21,*) XYp(p,1),XYp(p,2)
    end do
    close(21)

    ! ----------------------------------------
    filename  = 'maillage.'
    filename = trim(filename)//trim(cnum)

    open(21,file=filename)
    do j = 1,Mesh%nc
       do i = 1,Mesh%c_l(j)
          write(21,*) xn( Mesh%cell_list(j,i) ), yn( Mesh%cell_list(j,i) ), 0
       end do

       write(21,*)  xn( Mesh%cell_list(j,1) ), yn( Mesh%cell_list(j,1) ), 0
       write(21,*) ''

    end do
    close(21)

    ! ----------------------------------------
    filename  = 'lac.'
    filename = trim(filename)//trim(cnum)

    open(21,file=filename)
    do j = 1,1000
      call random_number(x)
      call random_number(y)
      x = 2 * x -1
      y = 2 * y -1
      xprime = x * sqrt(1 - (y**2)/2)
      yprime = y * sqrt(1 - (x**2)/2)
      write(21,*) (rayon- 2*1e-2) * xprime + x_lac, (rayon- 2*1e-2) * yprime + y_lac, 0
    end do
    close(21)

  end subroutine
end module our_module
