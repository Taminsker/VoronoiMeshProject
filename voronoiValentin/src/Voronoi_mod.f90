module Voronoi_mod

  use VPoint_mod
  use VEdge_mod
  use VEvent_mod
  use VParabola_mod

  implicit none
  save

  type(VPoint), dimension(:), pointer :: Vertices
  type(VEdge), dimension(:), pointer  :: Edges
