module mesh_data
  implicit none
  save
   integer, parameter :: d = selected_real_kind(P=15, R=50)

   type mesh_struct
     sequence

     ! cells
     integer                               :: numOfCells
     real(kind=d), dimension(:,:), pointer :: coorCent, coorMassC
     integer,      dimension(:),   pointer :: cellsList, areaOfCell

     ! nodes
     integer                               :: numOfGen, numOfNodes
     integer,      dimension(:,:), pointer :: genList, nodesList

     ! corners
     real(kind=d)                         :: xmax, xmin, ymax, ymin

   end type mesh_struct

end module mesh_data
