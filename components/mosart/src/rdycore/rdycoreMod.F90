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

  type(RDy)                               :: rdy_                      ! RDycore data structure

  PetscInt                                :: num_cells_owned           ! number of cells that locally owned
  PetscInt                                :: num_cells_global          ! total number of cells in the mesh
  PetscInt              , pointer         :: natural_id_cells_owned(:) ! natural IDs of cells that are locally owned

  integer               , public, pointer :: rdycore_pocn(:)           ! PE rank for each grid cell

  PetscReal             , pointer         :: total_runoff_data(:)      ! the water source to RDycore's SWE

  type(lnd2rdy_type)    , public          :: lnd2rdy_vars              ! data struture saving data sent from lnd to rdycore
  type(rdy_bounds_type) , public          :: rdy_bounds                ! bounds of grid cells

  integer               , public          :: iulog = 6

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
    PetscInt              :: g, myrank
    Vec                   :: owner_mpi, owner_seq
    PetscScalar, pointer  :: vec_ptr(:), rank(:)
    PetscInt   , pointer  :: int_ptr(:)
    IS                    :: is_from, is_to
    VecScatter            :: scatter
    character(len=100)    :: string
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

    allocate(natural_id_cells_owned(num_cells_owned))
    PetscCallA(RDyGetLocalCellNaturalIDs(rdy_, num_cells_owned, natural_id_cells_owned, ierr))

    PetscCallA(RDyGetNumGlobalCells(rdy_, num_cells_global, ierr))

    ! create a MPI and sequential Vec
    PetscCallA(VecCreateMPI(PETSC_COMM_WORLD, num_cells_owned, PETSC_DETERMINE, owner_mpi, ierr))
    PetscCallA(VecCreateSeq(MPI_COMM_SELF, num_cells_global, owner_seq, ierr))

    ! find the rank for the processor
    PetscCallA(mpi_comm_rank(PETSC_COMM_WORLD, myrank, ierr))

    allocate(rank(num_cells_owned))
    write(iulog,*)'myrank ',myrank
    rank(:) = myrank
    PetscCallA(VecSetValues(owner_mpi, num_cells_owned, natural_id_cells_owned, rank, INSERT_VALUES, ierr))
    deallocate(rank)
    PetscCallA(VecAssemblyBegin(owner_mpi, ierr))
    PetscCallA(VecAssemblyEnd(owner_mpi, ierr))

#if 0
       string = 'owner_mpi.txt'
       PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD, trim(string), viewer, ierr))
       PetscCallA(VecView(owner_mpi, viewer, ierr))
       PetscCallA(PetscViewerDestroy(viewer, ierr))
#endif

    ! create index set to scatter from the MPI Vec and to sequential Vec
    allocate(int_ptr(num_cells_global))
    do g = 1, num_cells_global
       int_ptr(g) = g - 1
    end do
    PetscCallA(ISCreateGeneral(PETSC_COMM_WORLD, num_cells_global, int_ptr, PETSC_COPY_VALUES, is_from, ierr))
    PetscCallA(ISCreateGeneral(PETSC_COMM_WORLD, num_cells_global, int_ptr, PETSC_COPY_VALUES, is_to, ierr))
    deallocate(int_ptr)

    ! create the VecScatter
    PetscCallA(VecScatterCreate(owner_mpi, is_from, owner_seq, is_to, scatter, ierr))
    PetscCallA(ISDestroy(is_from, ierr))
    PetscCallA(ISDestroy(is_to, ierr))

    ! scatter th data
    PetscCallA(VecScatterBegin(scatter, owner_mpi, owner_seq, INSERT_VALUES, SCATTER_FORWARD, ierr))
    PetscCallA(VecScatterEnd(scatter, owner_mpi, owner_seq, INSERT_VALUES, SCATTER_FORWARD, ierr))

#if 0
       write(string,*)myrank
       string = 'owner_seq_' // trim(adjustl(string)) // '.txt'
       PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_SELF, trim(string), viewer, ierr))
       PetscCallA(VecView(owner_seq, viewer, ierr))
       PetscCallA(PetscViewerDestroy(viewer, ierr))
#endif

    ! save PE number for each grid cell that MOSART will use
    allocate(rdycore_pocn(num_cells_global))

    ! determine the begin/end bounds of grid cells
    PetscCallA(VecGetArrayF90(owner_seq, vec_ptr, ierr))
    rdy_bounds%begg = 1
    do g = 1, num_cells_global
       rdycore_pocn(g) = int(vec_ptr(g))
       if (int(vec_ptr(g)) < myrank) then
          rdy_bounds%begg = rdy_bounds%begg + 1
       end if
    end do
    rdy_bounds%endg = rdy_bounds%begg - 1 + num_cells_owned
    PetscCallA(VecRestoreArrayF90(owner_seq, vec_ptr, ierr))

    ! allocate data structure for exchanging data from land to rdycore
    call lnd2rdy_vars%Init(rdy_bounds)

    ! free up memory
    PetscCallA(VecDestroy(owner_mpi, ierr))
    PetscCallA(VecDestroy(owner_seq, ierr))

    if (masterproc) then
       write(iulog,*)'RDycore model initialization completed'
    end if

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
    deallocate(rdycore_pocn)

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

