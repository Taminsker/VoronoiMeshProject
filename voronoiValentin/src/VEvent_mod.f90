module Vevent_mod
  use Vpoint_mod
  use Vparabola_mod

  implicit none
  save

  type Vevent
    sequence

    type(VPoint), pointer    :: point
    logical                  :: pe
    real*8                   :: y
    type(VParabola), pointer :: arch
  end type Vevent

contains
  subroutine VEvent_ctor(this, p, pev)
      type(VEvent), intent(inout)       :: this
      type(VPoint), pointer, intent(in) :: p
      logical, intent(inout)            :: pev

      this%point = p ; this%pe = pev ; this%y = p%y ; this%arch = 0;
  end subroutine VEvent_ctor

  function compareEvent(l, r)
    type(VEvent), pointer, intent(in) :: l
    type(VEvent), pointer, intent(in) :: r
    logical                           :: compareEvent

    compareEvent = (l%y < r%y)
  end function compareEvent

end module
