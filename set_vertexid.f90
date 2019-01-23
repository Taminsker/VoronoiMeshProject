subroutine set_vertexid(i,j,new_nid)
  !
  ! set vertex number to new_nid for all cells 
  !  that share vertex j of cell i
  !
  use mesh_data
  !
  implicit none
  !
  integer, intent(in) :: i,j,new_nid
  integer :: cell1, cell2,  vid1, vid2, nv
  !
  integer, dimension(100) :: iN, iNtest
  real(kind=real_acc), dimension(100) ::  xP,yP,xPtest,yPtest
  ! nPnt shold be the same as the dimension on the lines above
  integer, parameter :: nPnt=100
  !
  integer, external :: iGetNeibN
  !
  !----------------------------------------------------------
  !
  ! get cell info
  nv = iGetNeibN(i)
  call GetNbAr(i,iN,nPnt)
  call GetVAr(i,xP,yP,nPnt)
  !
  !
  Mesh%vertex_nid(i,j) = new_nid
  if(j.eq.1) Mesh%vertex_nid(i,nv+1) = new_nid
  !
  ! Now need to look for cells that also have this vertex
  !
  ! get cell id's of neighbour cells with this node
  cell1 = iN(j)
  if(cell1.gt.n) cell1 = 0
  cell2 = iN(j+1)
  if(cell2.gt.n) cell2 = 0
  ! if new node is boundary node then flag it
  if(cell1.eq.0.or.cell2.eq.0) Mesh%boundary(new_nid) = -1
  !
  if(cell1.gt.0) then
     ! for first cell identify vertex number
     ! know that cell vertices are held in counter-clockwise order
     nv = iGetNeibN(cell1)
     call GetNbAr(cell1,iNtest,nPnt)
     call GetVAr(cell1,xPtest,yPtest,nPnt)
     vid1=2
     do
        if(iNtest(vid1).eq.i) exit
        vid1 = vid1 + 1
     enddo
     vid1 = vid1 - 1
     Mesh%vertex_nid(cell1,vid1) = new_nid
     if(vid1.eq.1) Mesh%vertex_nid(cell1,nv+1) = new_nid
  endif
  !
  if(cell2.gt.0) then
     ! for second cell identify vertex number
     ! know that cell vertices are held in counter-clockwise order
     nv = iGetNeibN(cell2)
     call GetNbAr(cell2,iNtest,nPnt)
     call GetVAr(cell2,xPtest,yPtest,nPnt)
     vid2=1
     do
        if(iNtest(vid2).eq.i) exit
        vid2 = vid2 + 1
     enddo
     Mesh%vertex_nid(cell2,vid2) = new_nid
     if(vid2.eq.1) Mesh%vertex_nid(cell2,nv+1) = new_nid
  endif
  !
  return
  !
end subroutine set_vertexid

