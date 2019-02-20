module bibli_init
  use mesh_data
  implicit none

contains

  ! INITIALISATION OF PROBLEM TYPE
  !
  ! 0: random position of generators
  ! 1: uniform
  ! 2: -
  ! 3: -
  ! 4: Nothing
  ! 5: Annulus shape
  ! 6: full disk
  ! 7: -
  ! default: read a file of generators or mesh

  subroutine init_problem(problemNumber, mesh)
    implicit none

    !------------------------------------------------
    ! Parameters
    integer, intent(in) :: problemNumber
    type(mesh_struct), intent(inout) :: mesh

    real*8 :: randNumber

    !------------------------------------------------

    print*, '   ----> Enter in initialization problem #' , problemNumber

    select case(problemNumber)
    case(0)
      !----------------------
      ! Random
      integer :: cmpt = 0

      mesh%numOfGen = 250
      print*,'   How many generators ?'
      read*, mesh%numOfGen
      mesh%numOfGen = mesh%numOfGen + 4
      print*,' Add 4 generators at the four box corners'

      allocate(mesh%genList(1:mesh%numOfGen, 1:2))
      genList(1,1) = mesh%xmin ; genList(1,2) = mesh%ymin
      genList(2,1) = mesh%xmin ; genList(2,2) = mesh%ymax
      genList(3,1) = mesh%xmax ; genList(3,2) = mesh%ymin
      genList(4,1) = mesh%xmax ; genList(4,2) = mesh%ymax

      cmpt = cmpt + 4

      do i = 5,mesh%numOfGen
        call random_number(randNumber)
        genList(i,1) = mesh%xmin + abs(mesh%xmax - mesh%xmin) * randNumber
        call random_number(randNumber)
        genList(i,2) = mesh%ymin + abs(mesh%ymax - mesh%ymin) * randNumber
        cmpt = cmpt + 1
      end do
      if (.NOT. cmpt == mesh%numOfGen) then
        print*,'   ..Generators error numOfGen - cmpt = ',&
        & mesh%numOfGen- cmpt
      end if

      print*,'   ..Generators OK numOfGen',mesh%numOfGen

    case(1)
      !----------------------
      ! Uniform
      integer :: cmpt = 0

      mesh%numOfGen = 250
      print*,'   How many generators ?'
      read*, mesh%numOfGen
      mesh%numOfGen = mesh%numOfGen + 4

      real*8 :: dx = (xmx-xmn) / (sqrt(real(mesh%numOfGen)) )
      real*8 :: dy = (ymx-ymn) / (sqrt(real(mesh%numOfGen)) )
      real*8 :: ii = sqrt(real(mesh%numOfGen))
      print*,' Delta x, Delta y=',dx,dy,' Nc x Nc=',ii,'x',ii

      print*,' Add 4 generators at the four box corners'

      allocate(mesh%genList(1:mesh%numOfGen, 1:2))
      genList(1,1) = mesh%xmin ; genList(1,2) = mesh%ymin
      genList(2,1) = mesh%xmin ; genList(2,2) = mesh%ymax
      genList(3,1) = mesh%xmax ; genList(3,2) = mesh%ymin
      genList(4,1) = mesh%xmax ; genList(4,2) = mesh%ymax

      cmpt = cmpt + 4

      do q = 1,ii
        do p = 1,ii
          pp = pp + 1
          if( mod(q,2)==0 ) then
            XYp(pp,1) = xmn + dx/2.d0 + (p-1)*dx
            XYp(pp,2) = ymn + dy/2.d0 + (q-1)*dy
          else
            XYp(pp,1) = xmn + dx/2.d0 + (p-1)*dx -dx/5.d0
            XYp(pp,2) = ymn + dy/2.d0 + (q-1)*dy -dy/5.d0
          end if
          Np(pp) = 1
        end do
      end do
      p = pp
      npart = p
      ng = -1
      print*,'   ..Generators OK p=',p
    case(2)
      !
    case(3)
      !
    case(4)
      !
    case(5)
      !
    case(6)
      !
    case(7)
      !
    case default
      ! Read a file of generators
      filenamee = 'T.mesh'
      open(11,file=filenamee)
      read(11,*) npart, ng
      do i = 1,npart
        read(11,*) XYp(i,1),XYp(i,2)
      end do
      close(11)
      print*,'   ..Generators OK p=',npart
    end select

    ! Final computational bounds
    xmin = minval(XYp(:,1)); ymin = minval(XYp(:,2));
    xmax = maxval(XYp(:,1)); ymax = maxval(XYp(:,2));
    ! Final number of particles
    nnature = p
    !
    open(11,file='end_init.plt')
    do i = 1,npart
      write(11,*) XYp(i,1),XYp(i,2)
    end do
    close(11)

  end subroutine init_problem
