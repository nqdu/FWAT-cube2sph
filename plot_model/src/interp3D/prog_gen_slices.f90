program main 
  use mpi_f08
  use constants
  use interp3D_par
  use read_mesh3D
  use media_IO
  use locate_points 
  use interp3D 
  implicit none

  ! input filenames, external mesh path, model binary or dat path, points for interp, parameter name, outputfile name
  character(len=MAX_STR_LEN)  :: MESH_PATH,COOR_PATH,SUPRESS_UTM
  ! mpi init
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank,ier)

  ! get input args
  if(rank == 0) then 
    if(command_argument_count() /= 3) then 
      print *,'Usage :'
      stop 'need 3 parameters: mesh,coordinate, SUPRESS_UTM'
    else 
      call get_command_argument(1,MESH_PATH)
      call get_command_argument(2,COOR_PATH)
      MESH_PATH = trim(MESH_PATH)
      COOR_PATH = trim(COOR_PATH)
      call get_command_argument(3,SUPRESS_UTM)
      SUPRESS_UTM = trim(SUPRESS_UTM)
      write(*,'(a,a)') 'MESH_PATH = ',adjustl(trim(MESH_PATH))
      write(*,'(a,a)') 'COOR_PATH = ',adjustl(trim(COOR_PATH))    
      write(*,'(a,a)') 'SUPRESS_UTM= ',adjustl(trim(SUPRESS_UTM))   
    endif
  endif 
  ! bcast parameters to file
  call MPI_Bcast(MESH_PATH,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(COOR_PATH,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(SUPRESS_UTM,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  !setup gll points, weights, derivates
  
  !read external mesh
  call read_external_mesh_for_initialize(MESH_PATH)
  call read_external_mesh(MESH_PATH)
  !read points for interp
  call read_coor_get_lines(COOR_PATH)
  call read_coordinate(COOR_PATH,SUPRESS_UTM)

  ! locate coordinate (xi,eta,gammar) for each point
    call interp_initialize()
    !interp at (xi,eta,gamma) for each point
    !call interp()

    call write_profile(COOR_PATH,SUPRESS_UTM)

  
  call MPI_Barrier(MPI_COMM_WORLD)

  call MPI_FINALIZE(ier)
  
end program main 
