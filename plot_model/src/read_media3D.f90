program main 
  implicit none
  ! constants defined
  integer,PARAMETER   :: NGLLX = 5,NGLLY = NGLLX, NGLLZ = NGLLX 
  integer,PARAMETER   :: CUSTOM_REAL = 4,DOUBLE_REAL = 8

  ! filename
  character(len=256)  :: filename, paramfile,tempfile,LOCAL_PATH 
  
  ! useful arrays
  real(kind=CUSTOM_REAL),ALLOCATABLE   :: vp_read(:,:,:,:)  ! parameter
  integer,ALLOCATABLE :: ibool(:,:,:,:)   ! global index
  real(kind=CUSTOM_REAL),DIMENSION(:),ALLOCATABLE   :: x,y,z

  integer     :: nspec,nglob,i,j,k,ispec, rank,nproc,iglob,nspec_ire,ierr 
  integer,PARAMETER     :: IIN = 10,IOUT = 20

  ! get input args
  if(COMMAND_ARGUMENT_COUNT() /= 4) then 
    print *,'Usage :'
    stop './runthis LOCAL_PATH parameter nprocs out_name'
  else 
    call get_command_argument(3,paramfile)
    read(paramfile,*) nproc
    call get_command_argument(2,paramfile)
    tempfile = trim(paramfile)
    call get_command_argument(4,paramfile)
    !write(paramfile,'(a,a)')trim(paramfile),'.txt'
    call get_command_argument(1,LOCAL_PATH)
    LOCAL_PATH = trim(LOCAL_PATH)
    write(*,'(a,a)') 'LOCAL_PATH = ',adjustl(trim(LOCAL_PATH))
  endif

  ! open output file 
  open(IOUT,file=trim(paramfile))

  ! loop around each file
  do rank = 0,nproc-1 
    !read nspec and nglob
    write(filename,'(a,a,i6.6,a)') trim(LOCAL_PATH),'/proc',rank,'_external_mesh.bin'
    open(IIN,file=trim(filename),form='unformatted')
    read(IIN) nspec 
    read(IIN) nglob 
    read(IIN) nspec_ire
    
    ! allocate space
    ALLOCATE(vp_read(NGLLX,NGLLY,NGLLZ,nspec),ibool(NGLLX,NGLLY,NGLLZ,nspec),stat=ierr)
    if(ierr /=0) stop 'cannot allocate'
    ALLOCATE(x(nglob),y(nglob),z(nglob))

    ! read x,y,z,ibool
    x(:) = 0.0_CUSTOM_REAL; y(:) = 0.0_CUSTOM_REAL; z(:) = 0.0_CUSTOM_REAL
    ibool(:,:,:,:) = 0
    read(IIN) ibool
    read(IIN) x; read(IIN) y; read(IIN) z;
    close(IIN)

    print*,rank,nspec,nglob

    !read vp
    write(filename,'(a,a,i6.6,a,a,a)') trim(LOCAL_PATH),'/proc',rank,'_',trim(tempfile),'.bin'
    open(IIN,file=trim(filename),form='unformatted',iostat=ierr)
    if(ierr /=0 ) stop 'error readind vpfile '
    read(IIN) vp_read
    close(IIN)

    ! write to output file
    do ispec=1,nspec 
      do k=1,NGLLZ 
        do j=1,NGLLY 
          do i=1,NGLLX 
            iglob = ibool(i,j,k,ispec)
            utm_proj:block
              real(kind=DOUBLE_REAL)    :: rx,ry,rlon,rlat 
              rx = dble(x(iglob)); ry = dble(y(iglob))
              call utm_geo(rlon,rlat,rx,ry,10,1)
              write(IOUT,*) rlon,rlat,z(iglob) / 1000.,vp_read(i,j,k,ispec)
            end block utm_proj
          enddo
        enddo
      enddo
    enddo

    DEALLOCATE(x,y,z,ibool,vp_read)
  enddo

  ! close output file
  close(IOUT)
  
end program main 
