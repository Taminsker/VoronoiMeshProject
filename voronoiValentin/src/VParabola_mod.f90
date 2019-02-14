module VParabola_mod

  use VPoint_mod
  use VEdge_mod
  use VEvent_mod

  implicit none
  save

  type VParabola
    sequence
    logical                   :: isLeaf
    type(VPoint), pointer     :: site
    type(VEdge), pointer      :: edge
    type(VEvent), pointer     :: cEvent
    type(VParabola), pointer  :: parent
    type(VParabola), pointer  :: paraLeft
    type(VParabola), pointer  :: paraRight
  end type VParabola

contains
  subroutine VParabola_ctor(this)
    type(VParabola), intent(inout) :: this

    this%site = null()
    this%isLeaf = .FALSE.
    this%cEvent = null()
    this%edge = null()
    this%parent = null()
  end subroutine VParabola_ctor

  subroutine VParabola_ctor0(this, s)
    type(VParabola), intent(inout)        :: this
    type(VPoint), pointer, intent(intent) :: s

    this%site = s
    this%isLeaf = .TRUE.
    this%cEvent = null()
    this%edge = null()
    this%parent = null()
  end subroutine VParabola_ctor0

  subroutine SetLeft(this, p)
    type(VParabola), intent(inout)          :: this
    type(VParabola), intent(inout), pointer :: p

    this%paraLeft = p
    p%parent = this
  end subroutine SetLeft

  subroutine SetRight(this, p)
    type(VParabola), intent(inout)          :: this
    type(VParabola), intent(inout), pointer :: p

    this%paraRight = p
    p%parent = this
  end subroutine SetLeft

  function Left(this)
    type(VParabola), pointer, intent(inout) :: this
    type(VParabola), pointer                :: Left

    Left = this%paraLeft
  end function

  function Right(this)
    type(VParabola), pointer, intent(inout) :: this
    type(VParabola), pointer                :: Right

    Right = this%paraRight
  end function

  function GetLeft(p)
    type(VParabola), pointer, intent(inout) :: p
    type(VParabola), pointer                :: GetLeft

    GetLeft = GetLeftChild(GetLeftParent(p))
  end function

  function GetRight(p)
    type(VParabola), pointer, intent(inout) :: p
    type(VParabola), pointer                :: GetRight

    GetLeft = GetRightChild(GetRightParent(p))
  end function

  function GetLeftParent(p)
    type(VParabola), pointer, intent(inout) :: p
    type(VParabola), pointer                :: GetLeftParent

    type(VParabola), pointer                :: par = p%parent
    type(VParabola), pointer                :: pLast = p

    DO WHILE (Left(par) == pLast)
      if ( .NOT. par%parent ) then
        GetLeftParent = 0
        return
      end if
      pLast = par
      par = par%parent
    END DO

    GetLeftParent = par

  end function

  function GetRightParent(p)
    type(VParabola), pointer, intent(inout) :: p
    type(VParabola), pointer                :: GetRightParent

    type(VParabola), pointer                :: par = p%parent
    type(VParabola), pointer                :: pLast = p

    DO WHILE (Right(par) == pLast)
      if ( .NOT. par%parent ) then
        GetLeftParent = 0
        RETURN
      end if
      pLast = par
      par = par%parent
    END DO

    GetRightParent = par

  end function

  function GetLeftChild(p)
    type(VParabola), pointer, intent(inout) :: p
    type(VParabola), pointer                :: GetLeftChild

    if ( .NOT. p ) then
      GetLeftChild = 0;
      RETURN
    end if
    type(VParabola), pointer                :: par = Left(p)

    DO WHILE (.NOT. par%isLeaf)
      par = Right(par)
    END DO

    GetLeftChild = par

  end function

  function GetRightChild(p)
    type(VParabola), pointer, intent(inout) :: p
    type(VParabola), pointer                :: GetRightChild

    if ( .NOT. p ) then
      GetRightChild = 0;
      RETURN
    end if
    type(VParabola), pointer                :: par = Right(p)

    DO WHILE (.NOT. par%isLeaf)
      par = Left(par)
    END DO

    GetRightChild = par
  end function

end module VParabola_mod
