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
  character(len=MAX_STR_LEN)  :: MESH_PATH,MODEL_PATH,COOR_PATH,PAR_NAME,OUT_NAME, binORdat,SUPRESS_UTM,USE_INTERP
  ! mpi init
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank,ier)

  ! get input args
  if(rank == 0) then 
    if(command_argument_count() /= 8) then 
      print *,'Usage :'
      stop 'need 8 parameters: mesh,model,coordinate paths, parname, output filename, bin/dat SUPRESS_UTM USE_INTERP '
    else 
      call get_command_argument(1,MESH_PATH)
      call get_command_argument(2,MODEL_PATH)
      call get_command_argument(3,COOR_PATH)
      MESH_PATH = trim(MESH_PATH)
      COOR_PATH = trim(COOR_PATH)
      MODEL_PATH = trim(MODEL_PATH)
      call get_command_argument(4,PAR_NAME)
      call get_command_argument(5,OUT_NAME)
      call get_command_argument(6,binORdat)
      call get_command_argument(7,SUPRESS_UTM)
      call get_command_argument(8,USE_INTERP)
      PAR_NAME = trim(PAR_NAME)
      OUT_NAME = trim(OUT_NAME)
      binORdat = trim(binORdat) 
      write(*,'(a,a)') 'MESH_PATH = ',adjustl(trim(MESH_PATH))
      write(*,'(a,a)') 'COOR_PATH = ',adjustl(trim(COOR_PATH))
      write(*,'(a,a)') 'MODEL_PATH = ',adjustl(trim(MODEL_PATH))
      write(*,'(a,a)') 'MESH_PATH = ',adjustl(trim(PAR_NAME))
      write(*,'(a,a)') 'COOR_PATH = ',adjustl(trim(OUT_NAME))
      write(*,'(a,a)') 'MODEL_PATH = ',adjustl(trim(binORdat))     
      write(*,'(a,a)') 'SUPRESS_UTM= ',adjustl(trim(SUPRESS_UTM)) 
       write(*,'(a,a)') 'USE_INTERP = ',adjustl(trim(USE_INTERP))     
    endif
  endif 
  ! bcast parameters to file
  call MPI_Bcast(MESH_PATH,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(MODEL_PATH,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(COOR_PATH,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(PAR_NAME,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(OUT_NAME,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(binORdat,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(SUPRESS_UTM,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(USE_INTERP,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  !setup gll points, weights, derivates
  
  !read external mesh
  call read_external_mesh_for_initialize(MESH_PATH)
  call read_external_mesh(MESH_PATH)
  !read points for interp
  call read_coor_get_lines(COOR_PATH)
  call read_coordinate(COOR_PATH,SUPRESS_UTM)
  call read_binORdat_media(binORdat,MODEL_PATH,PAR_NAME)
  ! locate coordinate (xi,eta,gammar) for each point
  if(USE_INTERP=='true')then
    call interp_initialize()
    !interp at (xi,eta,gamma) for each point
    call interp()
    if(rank==0)call write_interped_media(OUT_NAME,SUPRESS_UTM)
    !for debug
  else
    call find_neareast()
    if(rank==0)call write_interped_media(OUT_NAME,SUPRESS_UTM)
  endif
  
  call MPI_Barrier(MPI_COMM_WORLD)

  call MPI_FINALIZE(ier)
  
end program main 
