module rdycoreMod

#include <petsc/finclude/petsc.h>

  use petsc
  use rdycore
  use RtmSpmd, only : mpicom_rof
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private

  type(RDy) :: rdy_

  public :: rdycore_init
  public :: rdycore_run
  public :: rdycore_final

contains

  !-----------------------------------------------------------------------
  subroutine rdycore_init()
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    character(len=1024)   :: config_file
    PetscErrorCode        :: ierr
 
    config_file = 'rdycore.yaml'
 
    ! set PETSc's communicator
    PETSC_COMM_WORLD = mpicom_rof
    PetscCallA(PetscInitialize(ierr))
 
    ! initialize subsystems
    PetscCallA(RDyInit(ierr))
 
    ! create rdycore and set it up with the given file
    PetscCallA(RDyCreate(PETSC_COMM_WORLD, config_file, rdy_, ierr))
    PetscCallA(RDySetup(rdy_, ierr))
 
  end subroutine rdycore_init

  !-----------------------------------------------------------------------
  subroutine rdycore_run()
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    use RtmTimeManager, only : get_step_size
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    real(r8)       :: dtime
    PetscErrorCode :: ierr

    dtime = get_step_size()

    ! Set the coupling time step
    PetscCallA(RDySetCouplingInterval(rdy_, dtime, ierr))

    ! Run the simulation to completion.
    PetscCallA(RDyAdvance(rdy_, ierr))
 
  end subroutine rdycore_run

  !-----------------------------------------------------------------------
  subroutine rdycore_final()
    !
    ! !DESCRIPTION:
    ! Destroy RDy object
    !
    ! !USES:
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    PetscErrorCode :: ierr

    PetscCallA(RDyDestroy(rdy_, ierr));
 
    PetscCallA(RDyFinalize(ierr));
 
  end subroutine rdycore_final
 
end module rdycoreMod

