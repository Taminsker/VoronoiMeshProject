module VEdge_mod
  use VPoint_mod

  implicit none
  save

  type VEdge
    sequence

    type(VPoint), pointer    :: start
    type(VPoint), pointer    :: end
    type(VPoint), pointer    :: direction
    type(VPoint), pointer    :: left
    type(VPoint), pointer    :: right
    real*8                   :: f
    real*8                   :: g
    type(VEdge), pointer    :: neighbour
  end type VEdge

contains
  subroutine VEdge_ctor(this, s, a, b)
    type(VEdge), intent(inout)            :: this
    type(VPoint), pointer, intent(inout)  :: s
    type(VPoint), pointer, intent(inout)  :: a
    type(VPoint), pointer, intent(inout)  :: b

    this%start = s
    this%left = a
    this%right = b
    this%neighbour = null()
    this%end = null()

    this%f = (b%x - a%x) / (a%y - b%y)
    this%g = s%y - this%f * s%x
    VPoint_ctor(this%direction, b%y - a%y, -(b%x - a%x))
  end subroutine
end module VEdge_mod
