module rdycoreMod

#include <petsc/finclude/petsc.h>

  use petsc
  use rdycore
  use RtmSpmd      , only : mpicom_rof, masterproc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
  use lnd2rdyType  , only : lnd2rdy_type
  use rdydecompMod , only : rdy_bounds_type

  implicit none

  private

  type(RDy)          :: rdy_

  PetscInt           :: num_cells_owned
  PetscReal, pointer :: total_runoff_data(:)

  type(lnd2rdy_type)   , public :: lnd2rdy_vars
  type(rdy_bounds_type), public :: rdy_bounds

  integer, public    :: iulog = 6

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

    call rdycore_setIO("rof_modelio.nml", iulog)
    if (masterproc) then
       write(iulog,*)'RDycore model initialization'
    end if

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
    allocate(total_runoff_data(num_cells_owned))

    if (masterproc) then
       write(iulog,*)'RDycore model initialization completed'
    end if

    rdy_bounds%begg = 1
    rdy_bounds%endg = num_cells_owned
    call lnd2rdy_vars%Init(rdy_bounds)

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
    PetscInt             :: t, nstep
    PetscReal            :: time_dn, time_up, cur_time, cur_rain
    PetscBool            :: found
    PetscInt             :: g, idx
    PetscErrorCode       :: ierr

    dtime    = get_step_size()
    nstep    = get_nstep()
    cur_time = (nstep-1)*dtime

    call get_curr_time_string(dateTimeString)
    if (masterproc) then
       write(*,*)'Beginning timestep of RDycore  : ',trim(dateTimeString)
       call shr_sys_flush(iulog)
    end if

    ! Set the water source term as the sum of surface and subsurface runoff
    do g = rdy_bounds%begg, rdy_bounds%endg
       idx = g - rdy_bounds%begg + 1
       total_runoff_data(idx) = lnd2rdy_vars%forc_qsur(g) + lnd2rdy_vars%forc_qsub(g)
    end do
    PetscCallA(RDySetWaterSourceForLocalCell(rdy_, num_cells_owned, total_runoff_data, ierr))

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

    ! close the logfile
    close(iulog)

    ! deallocate memory for rain data
    deallocate(total_runoff_data)

    ! destroy RDy object
    PetscCallA(RDyDestroy(rdy_, ierr));

    ! finalize
    PetscCallA(RDyFinalize(ierr));
 
  end subroutine rdycore_final
 
  !===============================================================================
  ! This subroutine is based on share/util/shr_file_mod.F90.
  ! The rof_modelio.nml file is opened and the value for the 'logfile' is read that
  ! would be something like the following:
  !
  ! logfile = "rof.log.230804-204952"
  !
  ! Based on the above-mentioned value, a logfile for RDycore is opend that would be
  ! "rdy.rof.log.230804-204952".
  !
  SUBROUTINE rdycore_setIO( nmlfile, funit)

    use shr_kind_mod
    use shr_sys_mod
    use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
    use shr_log_mod, only: s_loglev  => shr_log_Level
    use shr_log_mod, only: s_logunit => shr_log_Unit

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(len=*)    ,intent(in)  :: nmlfile  ! namelist filename
    integer(SHR_KIND_IN),intent(in)  :: funit    ! unit number for log file

    !EOP

    !--- local ---
    logical                :: exists   ! true if file exists
    character(SHR_KIND_CL) :: diri     ! directory to cd to
    character(SHR_KIND_CL) :: diro     ! directory to cd to
    character(SHR_KIND_CL) :: logfile  ! open unit 6 to this file
    character(SHR_KIND_CL) :: rdylogfile! open unit 6 to this file
    integer(SHR_KIND_IN)   :: unit     ! unit number
    integer(SHR_KIND_IN)   :: rcode    ! error code
    integer(SHR_KIND_IN)   :: l

    namelist / modelio / diri,diro,logfile

    !--- formats ---
    character(*),parameter :: subName = '(shr_file_setIO) '
    character(*),parameter :: F00   = "('(shr_file_setIO) ',4a)"
    character(*),parameter :: F01   = "('(shr_file_setIO) ',3a,i6)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !
    !-------------------------------------------------------------------------------

    diri = "."
    diro = "."
    logfile = ""

    inquire(file=nmlfile,exist=exists)

    if (.not. exists) then
       if (s_loglev > 0) write(s_logunit,F00) "file ",trim(nmlfile)," nonexistent"
       return
    else
       unit = shr_file_getUnit()
       open (unit,file=nmlfile,action="READ")
       read (unit,nml=modelio,iostat=rcode)
       close(unit)
       call shr_file_freeUnit( unit )
       if (rcode /= 0) then
          write(s_logunit,F01) 'ERROR: reading ',trim(nmlfile),': iostat=',rcode
          call shr_sys_abort(subName//" ERROR reading "//trim(nmlfile) )
       end if
    endif

    if (len_trim(logfile) > 0) then
       l = len(logfile)
       rdylogfile = logfile(5:l)
       open(funit,file=trim(diro)//"/rdy."//trim(rdylogfile))
    else
       if (s_loglev > 0) write(s_logunit,F00) "logfile not opened"
    endif

  END SUBROUTINE rdycore_setIO

end module rdycoreMod

