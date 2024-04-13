module rdydecompMod

  implicit none
  private

  public :: rdy_bounds_type

  type, public :: rdy_bounds_type
     integer :: begg, endg ! beginning and ending grid cell index
  end type rdy_bounds_type

end module rdydecompMod
