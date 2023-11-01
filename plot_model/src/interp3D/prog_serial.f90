program  main
  use constants
  use ieee_arithmetic,only : ieee_value,ieee_quiet_nan
  implicit none

  character(len=MAX_STR_LEN)  :: MODEL_PATH,COOR_PATH,PAR_NAME,OUT_NAME,temp
  integer :: myrank, i,j,k,nprocs,npts,ipts
  real,allocatable :: myprofile(:,:),param(:)
  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

  if(command_argument_count() /= 5) then 
    print *,'Usage :'
    stop 'need 5 parameters: MODEL,coordinate,PAR_NAME OUT_NAME nprocs'
  else 
    call get_command_argument(1,MODEL_PATH)
    call get_command_argument(2,COOR_PATH)
    MODEL_PATH = trim(MODEL_PATH)
    COOR_PATH = trim(COOR_PATH)
    call get_command_argument(3,PAR_NAME)
    call get_command_argument(4,OUT_NAME)
    PAR_NAME = trim(PAR_NAME)
    OUT_NAME = trim(OUT_NAME)
    write(*,'(a,a)') 'MODEL_PATH = ',adjustl(trim(MODEL_PATH))
    write(*,'(a,a)') 'COOR_PATH = ',adjustl(trim(COOR_PATH))    
    write(*,'(a,a)') 'PAR_NAME= ',adjustl(trim(PAR_NAME))
    write(*,'(a,a)') 'OUT_NAME= ',adjustl(trim(OUT_NAME))
    call get_command_argument(5,temp)
    read(temp,*) nprocs
    write(*,*) 'nprocs = ',nprocs
  endif

  ! read coordinates
  call read_coor_get_lines(COOR_PATH,npts)
  allocate(myprofile(9,npts),param(npts))
  param(:) = 0.
  open(10,file=trim(COOR_PATH))
  do i=1,npts 
    read(10,*)(myprofile(j,i),j=1,9)
    !if(i==10) print*,myprofile(:,i)
  enddo
  close(10)

  ! compute GLL points
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  print*,'serially interpolating ...'

  !$OMP parallel private(i,j,k,myrank,ipts)
  !$OMP parallel shared(MODEL_PATH,COOR_PATH,PAR_NAME,OUT_NAME,nprocs,npts,myprofile)
  !$OMP parallel shared(xigll,wxgll,yigll,wygll,zigll,wzgll,param)
  !$omp parallel do 
  do myrank = 0,nprocs-1
    tmp: block 
      character(len=MAX_STR_LEN) :: filename 
      integer :: ispec,nspec,ier
      double precision :: xi,eta,gamma 
      double precision :: hx(NGLLX),hpx(NGLLX),hy(NGLLY),hpy(NGLLY),hz(NGLLZ),hpz(NGLLZ)
      real(kind=CUSTOM_REAL),allocatable :: model_in(:,:,:,:)

      ! loop to find nspec for this rank
      nspec = 0
      do ipts = 1,npts 
        if(myrank == int(myprofile(7,ipts))) then 
          nspec = int(myprofile(9,ipts))
          exit 
        endif
      enddo
      !print*,myrank,nspec

      if(nspec == 0) cycle

      ! allocate space
      allocate(model_in(NGLLX,NGLLY,NGLLZ,nspec))

      ! read model
      write(filename,'(a,a,i6.6,a,a,a)') trim(MODEL_PATH),'/proc',myrank,'_',trim(PAR_NAME),'.bin'
      open(IIN,file=trim(filename),form='unformatted',iostat=ier)
      read(IIN) model_in
      close(IIN)

      ! interpolate
      do ipts = 1,npts 
        if(myrank /= int(myprofile(7,ipts))) cycle

        ! compute lagrange polynomial 
        xi = myprofile(4,ipts); eta = myprofile(5,ipts); gamma = myprofile(6,ipts)
        call lagrange_any(xi,NGLLX,xigll,hx,hpx)
        call lagrange_any(eta,NGLLY,yigll,hy,hpy)
        call lagrange_any(gamma,NGLLZ,zigll,hz,hpz)

        ispec = int(myprofile(8,ipts))
        do k=1,NGLLZ; do j=1,NGLLY; do i =1,NGLLX 
          param(ipts) = param(ipts) + hx(i) * hy(j) * hz(k) * model_in(i,j,k,ispec)
        enddo; enddo; enddo
      enddo

      deallocate(model_in)
    end block tmp  
  enddo
  !$OMP END PARALLEL DO

  ! save 
  where(param == 0.0)
    param = ieee_value(param,ieee_quiet_nan)
  endwhere
   
  open(IOUT,file=trim(OUT_NAME))
  do ipts=1,npts 
    write(IOUT,'(4(G0,1x))') myprofile(1:3,ipts), param(ipts)
  enddo
  close(IOUT)

end program  main

subroutine read_coor_get_lines(filename,npts)
  use constants
  implicit none

  character(len=MAX_STR_LEN) :: filename 
  ! local
  character(len=MAX_STR_LEN) :: line 
  integer :: npts,ier

  ! read until endoffile
  npts = 0
  open(IIN,file=filename)
  do 
    read(IIN,*,iostat=ier) line
    if(ier /= 0) exit 
    npts = npts + 1
  enddo
  close(IIN)
end subroutine read_coor_get_lines