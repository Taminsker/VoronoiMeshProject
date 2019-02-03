module mesh_data
  !
  !
  !
  ! This module contains all the mesh and connectivity variables required
  !  during a timestep, and functions used during the timestep
  !
  implicit none
  !
  save
  !
  ! Miminum accuracy for real numbers
  !
  integer,parameter :: real_acc = selected_real_kind(P=15,R=50)
  integer,parameter :: rr = real_acc
  integer,parameter :: d  = real_acc
  integer, parameter :: r_acc = real_acc
  integer, parameter :: maxe = 50
  integer, parameter :: maxn = 30000 * 3
  type Mesh_struct
     sequence
     ! Cell
     integer                                    :: nc
     real(kind=real_acc), dimension(:), pointer :: X_c, Y_c
     integer, dimension(:), pointer             :: c_l
     integer, dimension(:,:), pointer           :: cell_list
     integer, dimension(:,:), pointer           :: vertex_nid
     integer, dimension(:), pointer             :: cell_nid
     integer, dimension(:,:),pointer            :: neighbor
     ! Node
     integer                                    :: nn
     integer, dimension(:), pointer             :: n_l,boundary
     integer, dimension(:,:), pointer           :: node_list
     real(kind=real_acc), dimension(:), pointer :: X_n_n, Y_n_n
     integer, dimension(:), pointer             :: node_nid
     ! Corner
     integer                                             :: ncn
  end type Mesh_struct

  !----------------------------
  ! Mesh and Screen
  type(Mesh_struct)   :: Mesh
  !----------------------------
  !
  !----------------------------------------------------------
  !
  ! used when voronoi is invoked
  !
  ! Miminum accuracy for real numbers
  !
  !
  !
  integer, parameter :: MAXP=30000 * 3
  integer, parameter :: MAXTRI=2*MAXP+30
  integer, parameter :: MAXL=2000
  integer, parameter :: MAXG=MAXL
  integer, parameter :: MAXF=MAXL
  integer, parameter :: MAXDIR=6*MAXP+MAXF+30
  !
  !   +-----------+------------------------------------------------+
  !   |             |  Next triangle if we go over triangle ...      |
  !   | N of vert+-----------------------+------------------------+
  !   |             |  counter-clockwise    |       clockwise        |
  !   +-----------+-----------------------+------------------------+
  !   |    v[0]   |         t[1]          |          t[2]          |
  !   |    v[1]   |         t[2]          |          t[0]          |
  !   |    v[2]   |         t[0]          |          t[1]          |
  !   +-----------+-----------------------+------------------------+
  !
  type TriTable
     integer,dimension(0:3) :: v      ! Vertexes of triangle
     integer,dimension(0:3) :: t      ! References to neighb.triangle
  end type TriTable
  !
  type DirTable
     integer :: neib
     integer :: next
  end type DirTable
  !
  type(TriTable), dimension(0:MAXTRI) :: tri
  type(DirTable), dimension(0:MAXDIR) :: dNeib
  !
  real(kind=r_acc), dimension(0:MAXTRI) :: xv,yv
  real(kind=r_acc), dimension(0:MAXDIR) :: xVert,yVert
  real(kind=r_acc), dimension(0:MAXP+MAXF) :: x,y,xn,yn
  real(kind=r_acc), dimension(1:30000,1:3) :: FCT_ARRAY

  real(kind=r_acc) :: xmin,xmax,ymin,ymax
  real(kind=r_acc) :: eps = 100 * epsilon(eps)
  !
  integer :: n,nt,nFikt,dPusto,np,n1,n2,ng
  integer :: debug = 0, icycle
  integer, dimension(0:MAXTRI) :: ftri
  integer, dimension(0:MAXP+5) :: dRef
  integer, dimension(0:MAXF)   :: dType
  integer, dimension(0:MAXL)   :: l1,l2,l3
  integer, dimension(0:3)      :: LEFT,RIGHT
  !
  !-------------------------------------------------------


end	module mesh_data
