module Vpoint_mod
  implicit none
  save

  type VPoint
    sequence
    real*8 :: x
    real*8 :: y
  contains
    initial, pass :: VPoint_ctor
  end type VPoint

contains
  subroutine VPoint_ctor(this, vx, vy)
    type(VPoint), intent(inout) :: this
    real*8, intent(inout)       :: vx
    real*8, intent(inout)       :: vy

    this%x = vx ; this%y = vy
  end subroutine

end module
