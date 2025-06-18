module rdyconditionMod

#include <petsc/finclude/petsc.h>

  use petsc
  use rdycore
  use RtmSpmd      , only : mpicom_rof, masterproc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush

  implicit none
  private

#define DATASET_UNSET 0
#define DATASET_CONSTANT 1
#define DATASET_HOMOGENEOUS 2
#define DATASET_RASTER 3
#define DATASET_UNSTRUCTURED 4
#define DATASET_MULTI_HOMOGENEOUS 5

  type :: time_struct
    PetscInt :: year
    PetscInt :: month
    PetscInt :: day
    PetscInt :: hour
    PetscInt :: minute
  end type time_struct

  type :: UnstructuredDataset
    character(len=1024) :: dir
    character(len=1024) :: file
    character(len=1024) :: mesh_file
    character(len=1024) :: map_file

    PetscReal :: dtime_in_hour
    PetscInt  :: ndata_file

    Vec :: data_vec
    PetscScalar, pointer :: data_ptr(:)

    PetscInt :: ndata
    PetscInt :: stride

    type(time_struct)  :: start_date, current_date

    PetscInt :: header_offset

    PetscInt           :: mesh_nelements         ! number of cells or boundary edges in RDycore mesh
    PetscInt, pointer  :: data2mesh_idx(:)       ! for each RDycore element (cells or boundary edges), the index of the data in the unstructured dataset
    PetscReal, pointer :: data_xc(:), data_yc(:) ! x and y coordinates of data
    PetscReal, pointer :: mesh_xc(:), mesh_yc(:) ! x and y coordinates of RDycore elments

    PetscBool :: write_map_for_debugging         ! if true, write the mapping between the RDycore cells and the dataset for debugging
    PetscBool :: write_map                       ! if true, write the map between the RDycore cells and the dataset
    PetscBool :: read_map                        ! if true, read the map between the RDycore cells and the dataset

  end type UnstructuredDataset

  type, public :: BoundaryCondition
    PetscInt :: datatype

    type(UnstructuredDataset) :: unstructured

    PetscInt             :: ndata
    PetscInt             :: dirichlet_bc_idx
    PetscScalar, pointer :: data_for_rdycore(:)
  end type BoundaryCondition

  public :: CreateBoundaryConditionDataset
  public :: ParseBoundaryDataOptions
  public :: ApplyBoundaryCondition
contains

  ! Parse the command line options for the boundary condition dataset.
  ! The supported dataset types are:
  ! - spatially-homogeneous, temporally-varying
  ! - spatially and temporally varying unstructured dataset
  subroutine ParseBoundaryDataOptions(bc)
    !
    implicit none
    !
    type(BoundaryCondition) :: bc
    !
    PetscInt          :: dataset_type_count
    PetscBool         :: flg
    PetscBool         :: homogenous_bc_flag
    PetscBool         :: unstructured_bc_dir_flag
    PetscBool         :: unstructured_start_date_flag
    PetscInt, pointer :: date(:)
    PetscInt          :: ndate = 5
    PetscErrorCode    :: ierr

    bc%datatype                           = DATASET_UNSET
    bc%ndata                              = 0
    bc%dirichlet_bc_idx                   = -1

    dataset_type_count = 0

    ! parse information about unstructured boundary condition dataset
    PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-unstructured_bc_dir", bc%unstructured%dir, unstructured_bc_dir_flag, ierr))
    PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-unstructured_bc_write_map_for_debugging", bc%unstructured%write_map_for_debugging, flg, ierr))
    PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-unstructured_bc_write_map_file", bc%unstructured%map_file, bc%unstructured%write_map, ierr))
    PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-unstructured_bc_read_map_file", bc%unstructured%map_file, bc%unstructured%read_map, ierr))

    allocate(date(ndate))
    PetscCall(PetscOptionsGetIntArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-unstructured_bc_start_date", date, ndate, unstructured_start_date_flag, ierr))

    if (unstructured_start_date_flag) then
      dataset_type_count = dataset_type_count + 1
      if (ndate /= 5) then
        SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "-unstructured_bc_start_date should be in YY,MO,DD,HH,MM format")
      endif

      if (unstructured_bc_dir_flag .eqv. PETSC_FALSE) then
        SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "Need to specify path to unstructured BC data via -unstructured_bc_dir <dir>")
      endif

      PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-unstructured_bc_mesh_file", bc%unstructured%mesh_file, flg, ierr))
      if (flg .eqv. PETSC_FALSE) then
        SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "Need to specify the mesh file -unstructured_bc_mesh_file <file>")
      endif

      bc%datatype = DATASET_UNSTRUCTURED
      bc%unstructured%start_date%year   = date(1)
      bc%unstructured%start_date%month  = date(2)
      bc%unstructured%start_date%day    = date(3)
      bc%unstructured%start_date%hour   = date(4)
      bc%unstructured%start_date%minute = date(5)

      bc%unstructured%current_date%year   = date(1)
      bc%unstructured%current_date%month  = date(2)
      bc%unstructured%current_date%day    = date(3)
      bc%unstructured%current_date%hour   = date(4)
      bc%unstructured%current_date%minute = date(5)

    endif

    if (dataset_type_count > 1) then
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "More than one boundary condition type cannot be specified")
    endif

  end subroutine ParseBoundaryDataOptions

  subroutine FindDirichletBCID(rdy_, dirc_bc_idx, num_edges_dirc_bc, global_dirc_bc_idx, multiple_dirc_bcs_present)
    !
#include <finclude/rdycore.h>
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(RDy)      :: rdy_
    PetscInt       :: dirc_bc_idx, num_edges_dirc_bc
    PetscInt       :: global_dirc_bc_idx
    PetscBool      :: multiple_dirc_bcs_present
    !
    PetscInt       :: ibcond, nbconds, num_edges, bcond_type
    PetscErrorCode :: ierr

    dirc_bc_idx               = -1
    global_dirc_bc_idx        = -1
    num_edges_dirc_bc         = 0
    multiple_dirc_bcs_present = PETSC_FALSE

    PetscCallA(RDyGetNumBoundaryConditions(rdy_, nbconds, ierr))
    do ibcond = 1, nbconds
      PetscCallA(RDyGetNumBoundaryEdges(rdy_, ibcond, num_edges, ierr))
      PetscCallA(RDyGetBoundaryConditionFlowType(rdy_, ibcond, bcond_type, ierr))

      if (bcond_type == CONDITION_DIRICHLET) then
        if (dirc_bc_idx > -1) then
          multiple_dirc_bcs_present = PETSC_TRUE
        endif
        dirc_bc_idx       = ibcond
        num_edges_dirc_bc = num_edges
      endif
    enddo

    call MPI_Allreduce(dirc_bc_idx, global_dirc_bc_idx, 1, MPIU_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)

  end subroutine FindDirichletBCID

  ! Determines the filename of the dataset based on the directory and date
  ! The files are named as YYYY-MM-DD:HH-SS.<int32|int64>.bin
  subroutine DetermineDatasetFilename(dir, current_date, file)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    character(len=1024) :: dir
    type(time_struct)   :: current_date
    character(len=1024) :: file
    !
    PetscErrorCode      :: ierr

    write(file, '(A,"/",I4,"-",I2.2,"-",I2.2,":",I2.2,"-",I2.2,".",A,".bin")') trim(dir), current_date%year, current_date%month,  current_date%day, current_date%hour, current_date%minute, PETSC_ID_TYPE

  end subroutine DetermineDatasetFilename

  ! Opens the dataset file and reads the data into a PETSc Vec.
  ! The data is in the following format:
  !
  ! time_1 data_1
  ! time_2 data_2
  ! time_3 data_3
  subroutine opendata(filename, data_vec, ndata)
    implicit none
    character(*)   :: filename
    Vec            :: data_vec
    PetscInt       :: ndata

    PetscInt       :: size
    PetscViewer    :: viewer
    PetscErrorCode :: ierr

    PetscCallA(VecCreate(PETSC_COMM_SELF, data_vec, ierr))
    PetscCallA(PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, viewer, ierr))
    PetscCallA(VecLoad(data_vec, viewer, ierr));
    PetscCallA(PetscViewerDestroy(viewer, ierr));

    PetscCallA(VecGetSize(data_vec, size, ierr))
    ndata = size / 2

  end subroutine

  ! Opens the unstructured dataset file and reads the data into a PETSc Vec.
  subroutine OpenUnstructuredDataset(data, expected_data_stride)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(UnstructuredDataset) :: data
    PetscInt                  :: expected_data_stride
    !
    PetscInt                  :: size
    PetscErrorCode            :: ierr

    call DetermineDatasetFilename(data%dir, data%current_date, data%file)

    data%dtime_in_hour = 1.0
    data%ndata_file    = 1

    ! open the dataset file and read the data into a vector
    call opendata(data%file, data%data_vec, size)

    ! get the data pointer to the data in the vector
    PetscCallA(VecGetArray(data%data_vec, data%data_ptr, ierr))
    PetscCallA(VecGetSize(data%data_vec, size, ierr))

    data%header_offset = 2
    data%ndata         = data%data_ptr(1)
    data%stride        = data%data_ptr(2)

    if ((size - 2) / data%stride /= data%ndata) then
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "The number of data points in the unstructured dataset is not equal to the expected number of data points")
    endif

    if (data%stride /= expected_data_stride) then
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "The stride of the unstructured dataset is not equal to the expected stride")
    endif

  end subroutine OpenUnstructuredDataset

  ! Reads the x and y coordinates of the unstructured dataset mesh
  subroutine ReadUnstructuredDatasetCoordinates(data)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(UnstructuredDataset) :: data
    !
    PetscInt                 :: ndata, stride, offset, i
    PetscScalar, pointer     :: vec_ptr(:)
    Vec                      :: vec
    PetscViewer              :: viewer
    PetscErrorCode           :: ierr

    PetscCallA(VecCreate(PETSC_COMM_SELF, vec, ierr))
    PetscCallA(PetscViewerBinaryOpen(PETSC_COMM_SELF, data%mesh_file, FILE_MODE_READ, viewer, ierr))
    PetscCallA(VecLoad(vec, viewer, ierr))
    PetscCallA(PetscViewerDestroy(viewer, ierr))

    PetscCallA(VecGetArray(vec, vec_ptr, ierr))

    data%ndata = vec_ptr(1)
    allocate(data%data_xc(data%ndata))
    allocate(data%data_yc(data%ndata))

    stride = vec_ptr(2)
    if (stride /= 2) then
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "The stride of the unstructured dataset is not equal to 2")
    endif

    offset = 2;
    do i = 1, data%ndata
      data%data_xc(i) = vec_ptr(offset + (i - 1) * stride + 1)
      data%data_yc(i) = vec_ptr(offset + (i - 1) * stride + 2)
    enddo

    PetscCallA(VecRestoreArray(vec, vec_ptr, ierr))
    PetscCallA(VecDestroy(vec, ierr))

  end subroutine ReadUnstructuredDatasetCoordinates

  ! Creates the mapping between the RDycore cells and the unstructured dataset
  ! using nearest neighbor search
  subroutine CreateUnstructuredDatasetMap(data)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(UnstructuredDataset) :: data
    !
    PetscInt                 :: ndata, icell, i, count
    PetscReal                :: xc, yc, dx, dy, dist, min_dist

    ndata = data%ndata

    do icell = 1, data%mesh_nelements
      xc = data%mesh_xc(icell)
      yc = data%mesh_yc(icell)

      count = 0
      do i = 1, data%ndata
        count = count + 1
        dx = xc - data%data_xc(count)
        dy = yc - data%data_yc(count)
        dist = sqrt(dx * dx + dy * dy)
        if (i == 1) then
          min_dist = dist
          data%data2mesh_idx(icell) = count
        else
          if (dist < min_dist) then
            min_dist = dist
            data%data2mesh_idx(icell) = count
          endif
        endif
      enddo
    enddo

  end subroutine CreateUnstructuredDatasetMap

  ! Extracts the edge centroids from the RDycore mesh
  subroutine GetBoundaryEdgeCentroidsFromRDycoreMesh(rdy_, n, idx, xc, yc)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(RDy)               :: rdy_
    PetscInt                :: n, idx
    PetscReal, pointer      :: xc(:), yc(:)
    PetscErrorCode          :: ierr

    PetscCallA(RDyGetBoundaryEdgeXCentroids(rdy_, idx, n, xc, ierr))
    PetscCallA(RDyGetBoundaryEdgeYCentroids(rdy_, idx, n, yc, ierr))

  end subroutine GetBoundaryEdgeCentroidsFromRDycoreMesh

  ! After reading the unstructured dataset, this function
  ! postprocesses the data and sets up the mapping between
  ! the RDycore cells and the dataset.
  subroutine DoPostprocessForBoundaryUnstructuredDataset(rdy_, bc_dataset)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(RDy)               :: rdy_
    type(BoundaryCondition) :: bc_dataset
    !
    PetscInt               :: dirc_bc_idx, num_edges_dirc_bc
    PetscInt               :: global_dirc_bc_idx
    PetscBool              :: multiple_dirc_bcs_present

    call FindDirichletBCID(rdy_, dirc_bc_idx, num_edges_dirc_bc, global_dirc_bc_idx, multiple_dirc_bcs_present)

    ! do some sanity checking
    if (multiple_dirc_bcs_present) then
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "When BC file specified via -unstructured_bc_file argument, only one CONDITION_DIRICHLET can be present in the yaml")
    endif
    if (global_dirc_bc_idx == -1) then
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "No Dirichlet BC specified in the yaml file")
    endif

    bc_dataset%ndata            = num_edges_dirc_bc * 3
    bc_dataset%dirichlet_bc_idx = global_dirc_bc_idx

    allocate(bc_dataset%data_for_rdycore(bc_dataset%ndata))
    allocate(bc_dataset%unstructured%mesh_xc(num_edges_dirc_bc))
    allocate(bc_dataset%unstructured%mesh_yc(num_edges_dirc_bc))
    allocate(bc_dataset%unstructured%data2mesh_idx(bc_dataset%ndata))

    if (bc_dataset%ndata > 0) then
      bc_dataset%unstructured%mesh_nelements = num_edges_dirc_bc

      call GetBoundaryEdgeCentroidsFromRDycoreMesh(rdy_, num_edges_dirc_bc, global_dirc_bc_idx, bc_dataset%unstructured%mesh_xc, bc_dataset%unstructured%mesh_yc)

      call ReadUnstructuredDatasetCoordinates(bc_dataset%unstructured)

      ! set up the mapping between the dataset and boundary edges
      call CreateUnstructuredDatasetMap(bc_dataset%unstructured)

    endif

  end subroutine DoPostprocessForBoundaryUnstructuredDataset

  ! Set up boundary condition based on command line options
  subroutine CreateBoundaryConditionDataset(rdy_, bc_dataset)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(RDy)               :: rdy_
    type(BoundaryCondition) :: bc_dataset
    !
    PetscInt  :: expected_data_stride

    select case (bc_dataset%datatype)
    case (DATASET_UNSET)
      ! do nothing
    case (DATASET_UNSTRUCTURED)
      expected_data_stride = 3;
      call OpenUnstructuredDataset(bc_dataset%unstructured, expected_data_stride)
      call DoPostprocessForBoundaryUnstructuredDataset(rdy_, bc_dataset)
    case default
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "More than one boundary condition type cannot be specified")
    end select

  end subroutine CreateBoundaryConditionDataset

  ! Returns the number of days in a month
  subroutine day_of_month(year, month, ndays)
    PetscInt, intent(in) :: year, month
    PetscInt, intent(out) :: ndays

    select case (month)
      case (1, 3, 5, 7, 8, 10, 12)
        ndays = 31
      case (4, 6, 9, 11)
        ndays = 30
      case (2)
        if ((year / 400 == 0) .or. (year / 4 == 0) .and. (year / 100 /= 0)) then
          ndays = 29
        else
          ndays = 28
        end if
    end select

  end subroutine day_of_month

  ! Adds one hour to current date and updates the date if necessary
  subroutine IncrementDateByOneHour(current_date)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(time_struct) :: current_date
    !
    PetscInt          :: ndays

    current_date%hour = current_date%hour + 1

    if (current_date%hour == 24) then
      current_date%hour = 0
      current_date%day = current_date%day + 1

      call day_of_month(current_date%year, current_date%month, ndays)
      if (current_date%day > ndays ) then
        current_date%day = 1
        current_date%month = current_date%month + 1
        if (current_date%month > 12) then
          current_date%month = 1
          current_date%year = current_date%year + 1
        endif
      endif
    endif

  end subroutine IncrementDateByOneHour

  ! Closes the currently opened unstructured dataset and opens the next one
  subroutine OpenNextUnstructuredDataset(data)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(UnstructuredDataset) :: data
    !
    character(len=PETSC_MAX_PATH_LEN) :: outputString
    PetscInt                          :: tmpInt
    PetscErrorCode                    :: ierr

    PetscCallA(VecRestoreArray(data%data_vec, data%data_ptr, ierr))
    PetscCallA(VecDestroy(data%data_vec, ierr))

    call IncrementDateByOneHour(data%current_date)

    call OpenUnstructuredDataset(data, data%stride)

  end subroutine OpenNextUnstructuredDataset

  ! Sets spatially varying BC conditions from an unstructured dataset
  subroutine SetUnstructuredData(data, cur_time, num_values, data_for_rdycore)
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(UnstructuredDataset) :: data
    PetscReal                 :: cur_time
    PetscInt                  :: num_values
    PetscScalar, pointer      :: data_for_rdycore(:)
    !
    PetscInt                 :: ii, icell, idx, ndata_file, stride, offset

    if (cur_time / 3600.d0 >= (data%ndata_file) * data%dtime_in_hour) then
      ndata_file = data%ndata_file
      call OpenNextUnstructuredDataset(data)
      data%ndata_file = ndata_file + 1
    endif

    offset = data%header_offset;
    stride = data%stride

    do icell = 1, num_values
      idx = (data%data2mesh_idx(icell) - 1) * stride
      do ii = 1, stride
        data_for_rdycore((icell - 1) * stride + ii) = data%data_ptr(idx + ii + offset)
      enddo
    enddo

  end subroutine SetUnstructuredData

  ! Apply boundary condition to the RDycore object
  subroutine ApplyBoundaryCondition(rdy_, cur_time, bc_dataset)
#include <petsc/finclude/petsc.h>
#include <finclude/rdycore.h>
    !
    use rdycore
    use petsc
    !
    implicit none
    !
    type(RDy)               :: rdy_
    type(BoundaryCondition) :: bc_dataset
    PetscReal               :: cur_time
    !
    PetscInt, parameter     :: ndof = 3
    PetscErrorCode          :: ierr

    select case (bc_dataset%datatype)
    case (DATASET_UNSET)
      ! do nothing
    case (DATASET_UNSTRUCTURED)
      if (bc_dataset%ndata > 0) then
        call SetUnstructuredData(bc_dataset%unstructured, cur_time, bc_dataset%ndata / ndof, bc_dataset%data_for_rdycore)
        PetscCallA(RDySetFlowDirichletBoundaryValues(rdy_, bc_dataset%dirichlet_bc_idx, bc_dataset%ndata / ndof, ndof, bc_dataset%data_for_rdycore, ierr))
      endif
    case default
      SETERRA(PETSC_COMM_WORLD, PETSC_ERR_USER, "More than one boundary condition type cannot be specified")
    end select

  end subroutine ApplyBoundaryCondition

end module rdyconditionMod

