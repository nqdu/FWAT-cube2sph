module constants 

! constants defined
integer,PARAMETER   :: NGLLX = 5,NGLLY = NGLLX, NGLLZ = NGLLX 
integer,PARAMETER   :: CUSTOM_REAL = 4,DOUBLE_REAL = 8
integer,PARAMETER   :: UTM_PROJECTION_ZONE = 10
integer,parameter   :: MAX_STR_LEN = 256
integer,PARAMETER     :: IIN = 10,IOUT = 20

contains 

subroutine get_lines(filename,npts)
  implicit none
  
  character(len=MAX_STR_LEN) :: filename 
  integer,intent(inout)      :: npts 

  ! local
  character(len=MAX_STR_LEN) :: line 
  integer                    :: ierr 

  ! read until endoffile
  npts = 0
  open(IIN,file=filename)
  do 
    read(IIN,*,iostat=ierr) line
    if(ierr /= 0) exit 
    npts = npts + 1
  enddo
  close(IIN)

end subroutine get_lines

subroutine min_allprocs(vpmin_slice,distmin_slice,vpout,distout,npts)
  use mpi_f08
  implicit none
  integer,intent(in)   :: npts
  real(kind=CUSTOM_REAL),intent(in) :: vpmin_slice(npts),distmin_slice(npts)
  real(kind=CUSTOM_REAL),intent(inout) :: vpout(npts),distout(npts)

  ! local
  integer                 :: myrank,ierr,ipt
  real(kind=CUSTOM_REAL)  :: distin(2,npts),distall(2,npts),vmin(npts),dmin(npts)

  ! get current rank 
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank,ierr)

  ! find min distance in global
  distin(1,:) = distmin_slice(:)
  distin(2,:) = real(myrank,kind=CUSTOM_REAL)

  ! find minloc
  call MPI_ALLREDUCE(distin,distall,npts,MPI_2REAL,MPI_MINLOC,MPI_COMM_WORLD,ierr)
  vmin(:) = 0.0_CUSTOM_REAL
  dmin(:) = 0.0_CUSTOM_REAL

  do ipt=1,npts 
    if(myrank == int(distall(2,ipt))) then 
      vmin(ipt) = vpmin_slice(ipt)
      dmin(ipt) = distmin_slice(ipt)
    endif
  enddo
  
  call MPI_Reduce(vmin,vpout,npts,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(dmin,distout,npts,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

end subroutine min_allprocs

end module constants 


program main 
  use mpi_f08
  use constants
  implicit none

  ! filename
  character(len=MAX_STR_LEN)  :: LOCAL_PATH,INTFILE,parname,out_name,filename
  
  ! useful arrays
  real(kind=CUSTOM_REAL),ALLOCATABLE   :: vp_read(:,:,:,:)  ! parameter
  integer,ALLOCATABLE :: ibool(:,:,:,:)   ! global index
  real(kind=CUSTOM_REAL),DIMENSION(:),ALLOCATABLE   :: xstore,ystore,zstore

  ! interpolated coordinates
  real(kind=CUSTOM_REAL),ALLOCATABLE  :: xint(:),yint(:),zint(:),vpint(:)
  real(kind=CUSTOM_REAL),ALLOCATABLE  :: dist(:),distmin(:),vpmin(:)
  integer             :: npts 

  integer     :: nspec,nglob,i,j,k,ispec,iglob,nspec_ire,ierr,ipt 

  ! mpi parameters
  integer     :: rank,nproc

  ! temp
  real(kind=8)  :: rlon,rlat,rx,ry

  ! mpi init
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank,ierr)

  ! get input args
  if(rank == 0) then 
    if(command_argument_count() /= 4) then 
      print *,'Usage :'
      stop './runthis LOCAL_PATH INTFILE parname out_name'
      !call MPI_Abort(MPI_COMM_WORLD,1)
    else 
      call get_command_argument(1,LOCAL_PATH)
      LOCAL_PATH = trim(LOCAL_PATH)
      write(*,'(a,a)') 'LOCAL_PATH = ',adjustl(trim(LOCAL_PATH))

      call get_command_argument(2,INTFILE)
      INTFILE = trim(INTFILE)
      call get_command_argument(3,parname)
      call get_command_argument(4,out_name)      
    endif
  endif 
  
  ! bcast parameters to file
  call MPI_Bcast(LOCAL_PATH,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(INTFILE,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(parname,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)
  call MPI_Bcast(out_name,MAX_STR_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD)

  !read nspec and nglob
  write(filename,'(a,a,i6.6,a)') trim(LOCAL_PATH),'/proc',rank,'_external_mesh.bin'
  open(IIN,file=trim(filename),form='unformatted')
  read(IIN) nspec 
  read(IIN) nglob 
  read(IIN) nspec_ire
  
  ! allocate space
  ALLOCATE(vp_read(NGLLX,NGLLY,NGLLZ,nspec),ibool(NGLLX,NGLLY,NGLLZ,nspec),stat=ierr)
  if(ierr /=0) stop 'cannot allocate'
  ALLOCATE(xstore(nglob),ystore(nglob),zstore(nglob))

  ! read x,y,z,ibool
  read(IIN) ibool
  read(IIN) xstore; read(IIN) ystore; read(IIN) zstore;
  close(IIN)

  !read vp
  write(filename,'(a,a,i6.6,a,a,a)') trim(LOCAL_PATH),'/proc',rank,'_',trim(parname),'.bin'
  open(IIN,file=trim(filename),form='unformatted',iostat=ierr)
  if(ierr /=0 ) stop 'error readind vpfile '
  read(IIN) vp_read
  close(IIN)
  if(rank == 0)  then 
    write(filename,'(a,a,a,a,a)') trim(LOCAL_PATH),'/proc','_',trim(parname),'.bin'
    print*,'finished reading parameters '//trim(filename)
  endif

  ! read interpolated points
  call get_lines(INTFILE,npts)
  ALLOCATE(xint(npts),yint(npts),zint(npts),vpint(npts),&
            dist(npts),distmin(npts),vpmin(npts))
  open(IIN,file=INTFILE)
  do i=1,npts 
    read(IIN,*)rlon,rlat,zint(i)
    call utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,0) ! geo2utm
    xint(i) = real(rx,kind=CUSTOM_REAL)
    yint(i) = real(ry,kind=CUSTOM_REAL)
  enddo
  close(IIN)
  if(rank == 0)  then 
    print*,'finished reading coordinates'
  endif

  ! find the nearest point in this proc
  distmin(:) = huge(xint(1))
  do ispec = 1,nspec 
    do k=1,NGLLZ 
      do j=1,NGLLY 
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dist(:) = sqrt((xint(1:npts)-xstore(iglob))**2 &
                    +(yint(1:npts)-ystore(iglob))**2 &
                    +(zint(1:npts)-zstore(iglob))**2)
          do ipt=1,npts 
            if (dist(ipt) < distmin(ipt)) then 
              vpmin(ipt) = vp_read(i,j,k,ispec)
              distmin(ipt) = dist(ipt)
            endif 
          enddo
        enddo
      enddo
    enddo
  enddo
  if(rank == 0)  then 
    print*,'finished finding local distmin'
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! find the parameters for the nearest point in global
  call min_allprocs(vpmin,distmin,vpint,dist,npts)

  ! write file
  if(rank == 0) then 
    print *, 'Writing out gmt file ...'
    open(IOUT,file=out_name)
    do ipt=1,npts 
      rx = dble(xint(ipt)); ry = dble(yint(ipt))
      call utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,1) ! geo2utm
      !if(dist(ipt) > 1000.) vpint(ipt) = sqrt(-1.0)
      write(IOUT,'(5f12.5)')rlon,rlat,zint(ipt) * 0.001,vpint(ipt),dist(ipt)
    enddo
    close(IOUT)
  endif

  ! deallocate space
  DEALLOCATE(xstore,ystore,zstore,ibool,vp_read)
  DEALLOCATE(xint,yint,zint,vpint,dist,distmin,vpmin)

  ! close output file
  close(IOUT)

  call MPI_Barrier(MPI_COMM_WORLD)

  call MPI_FINALIZE(ierr)
  
end program main 
