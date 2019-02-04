module alloc

  use mesh_data

  implicit none

contains

  subroutine alloc_ini(Mesh)
    implicit none
    type(Mesh_struct),intent(inout) :: Mesh
    ! allocate memory for all ini values
    integer :: ierr
    !
    allocate (Mesh%boundary(maxn),stat=ierr)
    if(ierr.ne.0) then
       write(*,'("Allocation error in alloc - xn,yn,boundary")')
       
       stop
    endif
    !
    allocate (Mesh%c_l(maxn),Mesh%n_l(maxn),Mesh%cell_list(maxn,0:maxe),&
         &    Mesh%node_list(maxn,0:maxe), Mesh%neighbor(maxn,0:maxe),stat=ierr)
    if(ierr.ne.0) then
       write(*,'("Allocation error in alloc - c_l,n_l,cell_list,node_list")')
       stop
    endif
    !
    allocate (Mesh%vertex_nid(maxn,0:maxe),Mesh%cell_nid(maxn),Mesh%node_nid(maxn),stat=ierr)
    if(ierr.ne.0) then
       write(*,'("Allocation error in alloc - vertex_nid, cell_nid, node_nid")')
       stop
    endif
  end subroutine alloc_ini

    subroutine dealloc_ini_all(Mesh)
      implicit none
      type(Mesh_struct),intent(inout) :: Mesh
      ! deallocate memory for all ini values
      integer :: ierr
      deallocate ( Mesh%c_l,Mesh%n_l,Mesh%cell_list,Mesh%node_list,Mesh%neighbor,stat=ierr)
      deallocate ( Mesh%boundary,stat=ierr)
      deallocate ( Mesh%vertex_nid,Mesh%cell_nid,Mesh%node_nid,stat=ierr)
    end subroutine dealloc_ini_all

  end module alloc
