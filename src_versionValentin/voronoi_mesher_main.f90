program VORONOI

  use mesh_data
  use allocation
  use bibli_init
  use our_module

  implicit none
  !-----------------------------------------------------------
  ! Declarations
  !
  integer  :: problemNumber, maxIterationNumber, maxNumberGenerators
  type(mesh_data) :: mesh
  !-----------------------------------------------------------


  print*,'***************************************************************'
  print*,'*                                                             *'
  print*,'*   V O R O N O I   M E S H E R    2 D - New version          *'
  print*,'*                                                             *'
  print*,'***************************************************************'

 !-----------------------------------------------------------
 ! Problem choice
  print*,' Which initialisation of generators? (0:random,1:uniform,-1:from file)';
  read*, problemNumber

  !-----------------------------------------------------------
  ! Computational domain
  mesh%xmin = 0.0_d
  mesh%ymin = 0.0_d
  mesh%xmax = 1.0_d
  mesh%ymax = 1.0_d
  print*,' Domain to test: [', mesh%xmin, ',', mesh%xmax, &
          & '] x [', mesh%ymin, ',', mesh%ymax, ']'

  !-----------------------------------------------------------
  ! Define a maximal number of iteration
  maxIterationNumber = 10000

  !------------------------------------------------------------
  maxNumberGenerators = 50000
  print*,' Max number of generators : '
  print*,max_number_generators

  !------------------------------------------------------------
  call init_problem(problemNumber, mesh)


end program VORONOI
