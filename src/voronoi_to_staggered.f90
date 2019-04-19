subroutine voronoi_to_staggered(critical_length,critical_angle)
  !
  ! Convert output from Voronoi routines to staggered mesh
  !
  use mesh_data
  use alloc
  !
  implicit none
  !
  !
  integer :: i,j,k, counter,nn,nc
  integer :: new_nid,old_nid
  integer :: vid,cell
  integer :: nv, id
  real(kind=real_acc), intent(in) :: critical_length,critical_angle
  real(kind=real_acc) :: dist, check_dist, x_sum, y_sum, check
  real(kind=real_acc) :: dx1,dx2,dy1,dy2,len1,len2,angle
  !
  integer, dimension(100) :: iN
  real(kind=real_acc), dimension(100) ::  xP,yP
  ! nPnt should be the same as the dimension on the lines above
  integer, parameter :: nPnt=100
  !
  integer, external :: iGetNeibN
  !
  ! skipped_nid is a FIFO list of node id's that have been merged
  !  and so become unassigned
  logical :: skipped_node, not_found
  integer :: num_skipped_node, dim_skipped_nid
  integer, dimension(maxn) :: skipped_nid
  dim_skipped_nid=maxn
  !
  if(debug==1) print*,'    ** Subroutine voronoi_to_staggered => Generating staggered mesh from Voronoi cells'
  !
  ! Copy number of cell and node
  !
  Mesh%nn = nn
  Mesh%nc = nc
  !
  !
  ! initialise variables
  !
  skipped_node=.false.
  check_dist = critical_length**2
  if( check_dist > 0.0d0 ) print*,'CUT small edges at=',check_dist
  nn = 0
  num_skipped_node = 0
  do i=1,maxn
     Mesh%boundary(i) = 0
     do j=0,9
        Mesh%node_list(i,j) = 0
        Mesh%cell_list(i,j) = 0
        Mesh%vertex_nid(i,j) = 0
     enddo
  enddo
  !-------------------------------------------------------------------
  ! loop over cells to find node id for each cell vertex
  !
  if(debug==1) print*,'    First loop of voronoi_to_staggered n=',n,check_dist
  do i=1,n
     ! get number of vertices for cell n
     nv = iGetNeibN(i)
     ! get cell info
     call GetNbAr(i,iN,nPnt)
     call GetVAr(i,xP,yP,nPnt)
     ! loop over verticies of cell
     do j=1,nv
        ! check distance to next node in cells
        dist = (xP(j+1)-xP(j))**2 + (yP(j+1)-yP(j))**2
        if(dist.lt.check_dist) then
           !
           ! next vertex in cell is below critical distance
           !  give both verticies same node id
           !
           if(Mesh%vertex_nid(i,j).eq.0.and.Mesh%vertex_nid(i,j+1).eq.0) then
              ! neither vertex has already been assigned a node id
              if(skipped_node) then
                 new_nid = skipped_nid(num_skipped_node)
                 num_skipped_node = num_skipped_node - 1
                 if(num_skipped_node.eq.0) skipped_node = .false.
              else
                 nn = nn + 1
                 new_nid = nn
              endif
              ! assign vertices node number new_nid
              call set_vertexid(i,j,new_nid)
              if(j.lt.nv) then
                 call set_vertexid(i,j+1,new_nid)
              else
                 call set_vertexid(i,1,new_nid)
              endif
              !
              !
           else if(Mesh%vertex_nid(i,j).ne.0.and.Mesh%vertex_nid(i,j+1).eq.0) then
              ! current vertex already has node id, but next does not
              new_nid = Mesh%vertex_nid(i,j)
              call set_vertexid(i,j+1,new_nid)
              !
              !
           else if(Mesh%vertex_nid(i,j).eq.0.and.Mesh%vertex_nid(i,j+1).ne.0) then
              ! next vertex already has node id, but current does not
              new_nid = Mesh%vertex_nid(i,j+1)
              call set_vertexid(i,j,new_nid)
              !
              !
           else
              ! both vertices already have assigned node id
              !
              ! if both vertices already have same node id then do nothing
              if(Mesh%vertex_nid(i,j).ne.Mesh%vertex_nid(i,j+1)) then
                 ! assign smallest node id to both
                 skipped_node = .true.
                 num_skipped_node = num_skipped_node + 1
                 if(num_skipped_node.gt.dim_skipped_nid) then
                    write(*,*) 'Error, too many skipped nodes in subroutine mesh_gen (1)'
                    stop
                 endif
                 if(Mesh%vertex_nid(i,j).lt.Mesh%vertex_nid(i,j+1)) then
                    skipped_nid(num_skipped_node) = Mesh%vertex_nid(i,j+1)
                    Mesh%boundary(Mesh%vertex_nid(i,j+1)) = 0
                    new_nid = Mesh%vertex_nid(i,j)
                    if(j.lt.nv) then
                       call set_vertexid(i,j+1,new_nid)
                    else
                       call set_vertexid(i,1,new_nid)
                    endif
                 else
                    skipped_nid(num_skipped_node) = Mesh%vertex_nid(i,j)
                    old_nid = Mesh%vertex_nid(i,j)
                    new_nid = Mesh%vertex_nid(i,j+1)
                    ! Raph
                    Mesh%boundary(old_nid) = 0
                    call set_vertexid(i,j,new_nid,Mesh%vertex_nid)
                    ! PH, Raph
                    ! retrieve previous eventual skipped node
                    if( j==1 ) then
                       if(Mesh%vertex_nid(i,nv)==old_nid) then
                          skipped_nid(num_skipped_node) = Mesh%vertex_nid(i,nv)
                          new_nid = Mesh%vertex_nid(i,j+1)
                          ! Raph
                          Mesh%boundary(old_nid) = 0
                          call set_vertexid(i,nv,new_nid,Mesh%vertex_nid)
                       end if
                    else
                       if(Mesh%vertex_nid(i,j-1)==old_nid) then
                          skipped_nid(num_skipped_node) = Mesh%vertex_nid(i,j-1)
                          new_nid = Mesh%vertex_nid(i,j+1)
                          ! Raph
                          Mesh%boundary(old_nid) = 0
                          call set_vertexid(i,j-1,new_nid,Mesh%vertex_nid)
                       end if
                    end if
                 endif
              endif
           endif
           !
           !
        else
           ! next vertex in cell is above critical distance
           ! no merging
           if(Mesh%vertex_nid(i,j).eq.0) then
              if(skipped_node) then
                 new_nid = skipped_nid(num_skipped_node)
                 num_skipped_node = num_skipped_node - 1
                 if(num_skipped_node.eq.0) skipped_node = .false.
              else
                 nn = nn + 1
                 new_nid = nn
              endif
              ! assign vertices node number new_nid
              call set_vertexid(i,j,new_nid)
           endif
        endif
     enddo
     ! add first node to end of list
     Mesh%vertex_nid(i,nv+1) = Mesh%vertex_nid(i,1)
  enddo
  !
  ! finished loop over cells to find node id for each cell vertex
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  ! Now generate lists
  ! cell_list: for each node the cells that it is attached to
  ! node_list: for each cell the nodes that difine its boundary
  ! edge_list: for each cell the adjacent cells
  !
  do i=1,nn
     Mesh%n_l(i) = 0
  enddo
  !
  if(debug==1) print*,'      Phase 1 - done'
  if(debug==1) print*,'      Nc=',n
  !
  nc = n
  do i=1,nc
     Mesh%c_l(i) = 0
     !e_l(i) = 0
     j=1
     !
     Mesh%c_l(i)                    = Mesh%c_l(i) +1
     Mesh%cell_list(i,Mesh%c_l(i))  = Mesh%vertex_nid(i,j)
     Mesh%n_l(Mesh%vertex_nid(i,j)) = Mesh%n_l(Mesh%vertex_nid(i,j)) +1
     Mesh%node_list( Mesh%vertex_nid(i,j),Mesh%n_l(Mesh%vertex_nid(i,j)) ) = i
     !
     nv = iGetNeibN(i)
     call GetNbAr(i,iN,nPnt)
     do
        j=j+1
        if(j.gt.nv) exit
        !
        if(Mesh%vertex_nid(i,j).ne.Mesh%vertex_nid(i,j-1).and.Mesh%vertex_nid(i,j).ne.Mesh%vertex_nid(i,1)) then
           Mesh%c_l(i)                   = Mesh%c_l(i) +1
           Mesh%cell_list(i,Mesh%c_l(i)) = Mesh%vertex_nid(i,j)
           !
           Mesh%n_l(Mesh%vertex_nid(i,j))= Mesh%n_l(Mesh%vertex_nid(i,j)) + 1
           Mesh%node_list(Mesh%vertex_nid(i,j),Mesh%n_l(Mesh%vertex_nid(i,j))) = i
        endif
     enddo
     Mesh%cell_list(i,0)             = Mesh%cell_list(i,Mesh%c_l(i))
     Mesh%cell_list(i,Mesh%c_l(i)+1) = Mesh%cell_list(i,1)
     ! get edge lists
     !do j=1,nv
     !if(Mesh%vertex_nid(i,j).ne.Mesh%vertex_nid(i,j+1)) then
     !e_l(i) = e_l(i) + 1
     !edge_list(i,e_l(i)) = iN(j+1)
     !endif
     !enddo
  enddo
  ! finished generating cell, node and edge lists
  !------------------------------------------------------------
  !------------------------------------------------------------
  ! Now calculate nodal coordinates
  !
  if(debug==1) print*,'      Phase 2 - list generation - done'
  if(debug==1) print*,'      Nn=',nn
  do i=1,nn
     ! prevent division by zero when considering skipped node
     not_found=.true.
     if(skipped_node) then
        do j=1,num_skipped_node
           if(skipped_nid(j).eq.i) not_found=.false.
        enddo
     endif
     !
     if(not_found) then
        counter = 0
        x_sum = 0.0_d
        y_sum = 0.0_d
        ! if it is a boundary node then move internal verticies to boundary
        if(Mesh%boundary(i)/=0) then
           ! identify boundary vertex
           do j=1,Mesh%n_l(i)
              nv = iGetNeibN(Mesh%node_list(i,j))
              !id = Mesh%node_list(i,j)
              call GetNbAr(Mesh%node_list(i,j),iN,nPnt)
              call GetVAr(Mesh%node_list(i,j),xP,yP,nPnt)
              do k=1,nv
                 if(Mesh%vertex_nid(Mesh%node_list(i,j),k).eq.i) then
                    if(iN(k).eq.0.or.iN(k+1).eq.0) then
                       x_sum = xP(k)
                       y_sum = yP(k)
                    endif
                    ! New stuff added by PH, Misha and Raph June 09
                    if(iN(k)>nc.or.iN(k+1)>nc) then
                       x_sum = xP(k)
                       y_sum = yP(k)
                    endif
                 endif
              enddo
           enddo
           xn(i) = x_sum
           yn(i) = y_sum
        else
           ! otherwise calculate average coordinate
           do j=1,Mesh%n_l(i)
              nv = iGetNeibN(Mesh%node_list(i,j))
              call GetVAr(Mesh%node_list(i,j),xP,yP,nPnt)
              do k=1,nv
                 if(Mesh%vertex_nid(Mesh%node_list(i,j),k).eq.i) then
                    counter = counter + 1
                    x_sum = x_sum + xP(k)
                    y_sum = y_sum + yP(k)
                 endif
              enddo
           enddo
           xn(i) = x_sum/counter
           yn(i) = y_sum/counter
        endif
     endif
  end do
  ! finished calculating nodal coordinates
  !-----------------------------------------------------------
  !open(24,file='mesh.plt')
  !do cell = 1,nc
  !   do i = 1,Mesh%c_l(cell)
  !      write(24,*) xn( Mesh%cell_list(cell,i) ),yn( Mesh%cell_list(cell,i) )
  !   end do
  !   i=1
  !   write(24,*) xn( Mesh%cell_list(cell,i) ),yn( Mesh%cell_list(cell,i) )
  !   write(24,*) ''
  !end do
  !close(24)

  if(debug==1) print*,'      Phase 3 - Node list - done'
  !-----------------------------------------------------------
  ! Search for redundant nodes
  !  these are boundary nodes that are only associated with one cell
  !  where the adjoining edges meet at close to 180 degrees
  check = cos(critical_angle)
  do i=1,nn
     if(Mesh%n_l(i).eq.1) then
        cell = Mesh%node_list(i,1)
        vid = 0
        ! check angle of edges
        do j=1,Mesh%c_l(cell)
           if(Mesh%cell_list(cell,j).eq.i) then
              vid = j
              exit
           endif
        enddo
        dx1 = xn(i) - xn(Mesh%cell_list(cell,vid-1))
        dy1 = yn(i) - yn(Mesh%cell_list(cell,vid-1))
        dx2 = xn(Mesh%cell_list(cell,vid+1)) - xn(i)
        dy2 = yn(Mesh%cell_list(cell,vid+1)) - yn(i)
        len1 = sqrt(dx1*dx1 + dy1*dy1)
        len2 = sqrt(dx2*dx2 + dy2*dy2)
        angle = (dx1*dx2 + dy1*dy2)/(len1*len2)
        if(angle.gt.check) then
           if(debug==1) print*,' Small angle'
           ! small angle, get rid of node
           skipped_node = .true.
           num_skipped_node = num_skipped_node + 1
           if(num_skipped_node.gt.dim_skipped_nid) then
              write(*,*) 'Error, too many skipped nodes in subroutine mesh_gen (2)'
              stop
           endif
           skipped_nid(num_skipped_node) = i
           ! remove from cell list
           do j=1,Mesh%c_l(cell)
              if(j.gt.vid) then
                 Mesh%cell_list(cell,j-1) = Mesh%cell_list(cell,j)
                 !edge_list(cell,j-1) = edge_list(cell,j)
              endif
           enddo
           Mesh%c_l(cell) = Mesh%c_l(cell) - 1
           Mesh%cell_list(cell,Mesh%c_l(cell)+1) = Mesh%cell_list(cell,1)
           !e_l(cell) = e_l(cell) - 1
        endif
     endif
  enddo
  !finished search for redundant nodes
  !---------------------------------------------------------
  !do i=1,nc
  ! do j=1,Mesh%c_l(i)
  !  write(46,*) i,j,Mesh%cell_list(i,j)
  ! enddo
  !enddo
  !---------------------------------------------------------
  ! Sort lists so that there are no skipped node id's
  !
  if(debug==1) print*,'      Phase 4 - done'
  if(debug==1) print*,'      NN=',nc,nn,' Skipnode=',skipped_node,num_skipped_node
  !
  if(skipped_node) then
     counter = 0
     do i=1,nn
        not_found = .true.
        do j=1,num_skipped_node
           if(skipped_nid(j).eq.i) not_found=.false.
        enddo
        if(not_found) then
           counter = counter + 1
           Mesh%node_nid(counter) = i
           Mesh%cell_nid(i) = counter
        else
           Mesh%cell_nid(i) = 0
        endif
     enddo
     if(debug==1) print*,'      Phase 5 - Skip node 1/2 - done'
     !
     ! move node variables
     nn = nn - num_skipped_node
     do i=1,nn
        xn(i) = xn(Mesh%node_nid(i))
        yn(i) = yn(Mesh%node_nid(i))
        Mesh%n_l(i) = Mesh%n_l(Mesh%node_nid(i))
        Mesh%boundary(i) = Mesh%boundary(Mesh%node_nid(i))
        do j=1,Mesh%n_l(i)
           Mesh%node_list(i,j) = Mesh%node_list(Mesh%node_nid(i),j)
        enddo
     enddo
     ! revised cell list
     do i=1,nc
        do j=0,Mesh%c_l(i)+1
           Mesh%cell_list(i,j) = Mesh%cell_nid(Mesh%cell_list(i,j))
        enddo
     enddo
     if(debug==1) print*,'      Phase 5 - Skip node 2/2 - done'
  endif
  !finished sorting lists
  !---------------------------------------------------------

  ! Final number of points/cells
  Mesh%nn = nn
  Mesh%nc = nc

  ! SCREEN PRINTS (DEBUG)
  !      write(*,*) '1<= i <= Mesh%nn, 1<= j <= Mesh%c_l(i) , Mesh%cell_list(i,j)'
  !      do i=1,nc
  !       !do j=1,Mesh%c_l(i)
  !        write(*,*) 'Cell ',i,Mesh%c_l(i),' nodes asso:',Mesh%cell_list(i,:Mesh%c_l(i)+1)
  !       !enddo
  !      enddo
  !      write(*,*) '1<= i <= Mesh%nc, 1<= j <= Mesh%n_l(i) , Mesh%node_list(i,j)'
  !      do i=1,nn
  !       !do j=1,Mesh%n_l(i)
  !        write(*,*) 'Node ',i,Mesh%n_l(i),' cells asso:',Mesh%node_list(i,:Mesh%n_l(i)+1),' Xi,Yi=',xn(i),yn(i)
  !       !enddo
  !      enddo

  if(debug==1) print*,'    ** End Subroutine voronoi_to_staggered ncells,nnodes=',Mesh%nc,Mesh%nn

end subroutine voronoi_to_staggered
