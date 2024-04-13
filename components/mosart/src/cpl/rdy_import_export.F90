module rdy_import_export

  use shr_kind_mod , only : r8 => shr_kind_r8
  use rdydecompMod , only : rdy_bounds_type
  use lnd2rdyType  , only : lnd2rdy_type
  use rdy_cpl_indices

  implicit none
  private

  public :: rdy_import_mct

contains

  !------------------------------------------------------------------------
  subroutine rdy_import_mct(rdy_bounds, x2rd, lnd2rdy_vars)
    !
    ! !DESCRIPTION:
    ! Imports data from the coupler for the RDycore
    !
    implicit none
    !
    ! !ARGUMENTS:
    type(rdy_bounds_type) , intent (in)    :: rdy_bounds
    real(r8)              , intent (in)    :: x2rd(:,:)
    type(lnd2rdy_type)    , intent (inout) :: lnd2rdy_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g, idx                         ! indices
    real(r8) :: runoff_unit_conversion = 1.d-3 ! [mm/s] to [m/s]

    do g = rdy_bounds%begg, rdy_bounds%endg

       idx = g - rdy_bounds%begg + 1

       lnd2rdy_vars%forc_qsur(g) = x2rd(index_x2rdy_Flrl_qsur, idx) * runoff_unit_conversion
       lnd2rdy_vars%forc_qsub(g) = x2rd(index_x2rdy_Flrl_qsub, idx) * runoff_unit_conversion

    end do

  end subroutine rdy_import_mct

end module rdy_import_export
