!*****************************************************
program mpicalldynamo
  use mpi
  use dynamo
  use input_params
  use global_input_parameters
!
  implicit none 
!
  logical, parameter :: master_participate= .false. 
  integer :: igal, jgal, nmygals, nslaves, order, flag
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, parameter :: master=0
  integer, allocatable, dimension(:) :: mygals
  character(len=32) :: command_argument

!*****************************************************
  real :: tstart,tfinish
  integer :: rank, nproc, ierr, rc, len
  character*(MPI_MAX_PROCESSOR_NAME) name
! 
  call MPI_INIT(ierr)
  if (ierr/= MPI_SUCCESS) then
    print*,'Error starting MPI program. Terminating.'
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) !Get the rank of the processor this thread is running on
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr) !Get the number of processors this job
  call MPI_GET_PROCESSOR_NAME(name, len, ierr) !Get the name of this processor (usually the hostname)
  if (ierr/= MPI_SUCCESS) then
    print*,'Error getting processor name. Terminating.'
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif
  
  call get_command_argument(1, command_argument)
  if (len_trim(command_argument) == 0) then
    write(*,*) 'No parameter file provided. Using standard global parameters.'
  else
    call read_global_parameters(trim(command_argument))
  endif
  
  !*****************************************************
  !
  ! Assign a chunk of galaxies to each processor
  !
  if (rank== 0) then !Only the master (rank 0)
    tstart= MPI_wtime()
  endif
  !
  if (master_participate) then
    nslaves=nproc
    order=rank
  else
    nslaves=nproc-1
    order=rank-1
  endif
  !
  nmygals=0
  do igal=1,ngals
    if (mod(igal,nslaves) == order) then
      nmygals=nmygals+1
    endif
  enddo
  allocate(mygals(nmygals))
  jgal=1
  do igal=1,ngals
    if (mod(igal,nslaves) == order) then
      mygals(jgal)=igal
      jgal=jgal+1
    endif
  enddo

  print*,'rank=',rank,'    mygals=',mygals !processor id, list of galaxies for that processor

  ! Call dynamo code for each galaxy in mygals
  if (nmygals > 0) then
    do igal=1,nmygals
      call dynamo_run(info, igal, flag)
      print*,'flag',flag,'obtained by processor',rank,'for galaxy',mygals(igal)
    enddo
  endif

  call MPI_FINALIZE(ierr) !Tell the MPI library to release all resources it is using
  
  if (rank == 0) then !Only the master (rank 0)
    tfinish= MPI_wtime()
    print*,'total wall time in seconds=',tfinish-tstart
  endif
end program mpicalldynamo
!*****************************************************
