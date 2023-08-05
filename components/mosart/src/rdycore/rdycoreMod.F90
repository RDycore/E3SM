module rdycoreMod

#include <petsc/finclude/petsc.h>

  use petsc
  use rdycore
  use RtmSpmd      , only : mpicom_rof, masterproc
  use RtmVar       , only : iulog
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush

  implicit none

  private

  type(RDy) :: rdy_
  Vec       :: rain_timeseries
  PetscInt  :: rain_stride
  PetscInt  :: rain_nstep
  PetscInt  :: num_cells_owned
  PetscReal, pointer :: rain_data(:)

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
    PetscViewer           :: viewer
    PetscInt              :: size
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
 
    ! allocate memory for grid-level rain data
    PetscCallA(RDyGetNumLocalCells(rdy_, num_cells_owned, ierr))
    allocate(rain_data(num_cells_owned))

    ! Read the rain vector that has the following format:
    !
    !  time_1 rain_value_1
    !  time_2 rain_value_2
    !
    PetscCallA(VecCreate(PETSC_COMM_WORLD, rain_timeseries, ierr))
    PetscCallA(PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'rain.bin', FILE_MODE_READ, viewer, ierr))
    PetscCallA(VecLoad(rain_timeseries, viewer, ierr))
    PetscCallA(PetscViewerDestroy(viewer, ierr))

    ! Determine the number of rain
    rain_stride = 2
    PetscCallA(VecGetSize(rain_timeseries, size, ierr))
    rain_nstep = size/rain_stride

  end subroutine rdycore_init

  !-----------------------------------------------------------------------
  subroutine rdycore_run()
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    use RtmTimeManager, only : get_step_size, get_nstep, get_curr_time_string
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    character(len=256)   :: dateTimeString
    real(r8)             :: dtime
    PetscScalar, pointer :: rain_p(:)
    PetscInt             :: t, nstep
    PetscReal            :: time_dn, time_up, cur_time, cur_rain
    PetscBool            :: found
    PetscErrorCode       :: ierr

    dtime    = get_step_size()
    nstep    = get_nstep()
    cur_time = (nstep-1)*dtime

    call get_curr_time_string(dateTimeString)
    if (masterproc) then
       write(*,*)'Beginning timestep of RDycore  : ',trim(dateTimeString)
       call shr_sys_flush(iulog)
    end if


    ! Find the current rainfall
    PetscCallA(VecGetArrayF90(rain_timeseries, rain_p, ierr))
    found = PETSC_FALSE

    do t = 1, rain_nstep-1
       time_dn = rain_p((t-1)*rain_stride + 1)
       time_up = rain_p((t-1)*rain_stride + 3)
       if (cur_time >= time_dn .and. cur_time < time_up) then
          found = PETSC_TRUE
          cur_rain = rain_p((t-1)*rain_stride + 2)
          exit
       end if
    end do

    if (.not.found) then
       cur_rain = rain_p((nstep-1)*rain_stride + 2)
    end if

    PetscCallA(VecRestoreArrayF90(rain_timeseries, rain_p, ierr))

    ! Set spatially homogeneous rainfall for all grid cells
    rain_data(:) = cur_rain
    PetscCallA(RDySetWaterSource(rdy_, rain_data, ierr))

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

    ! deallocate memory for rain data
    deallocate(rain_data)
    PetscCallA(VecDestroy(rain_timeseries, ierr))

    ! destroy RDy object
    PetscCallA(RDyDestroy(rdy_, ierr));

    ! finalize
    PetscCallA(RDyFinalize(ierr));
 
  end subroutine rdycore_final
 
end module rdycoreMod

