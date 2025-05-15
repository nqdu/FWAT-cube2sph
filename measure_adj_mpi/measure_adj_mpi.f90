subroutine read_par_file_mpi(myrank,fstart0,fend0,tt,dtt,nn,chan)

  use ma_constants
  use ma_variables
  use ma_sub2
  use ascii_rw
  implicit none
  integer,intent(in) :: myrank 
  double precision, intent(out) :: fstart0,fend0,tt,dtt
  integer, intent(out) :: nn
  character(len=10), intent(out) :: chan
  integer :: ios

  ! input file MEASUREMENT.PAR -- see write_par_file.pl for details

  OUT_DIR = 'OUTPUT_FILES'   ! default

  open(10,file='MEASUREMENT.PAR',status='old',iostat=ios)
  if ( ios /= 0) stop 'Error opening MEASUREMENT.PAR file'
  read(10,*) tt,dtt,nn
  read(10,*) imeas0
  read(10,*) chan
  read(10,*) TLONG, TSHORT
  read(10,*) RUN_BANDPASS
  read(10,*) DISPLAY_DETAILS
  read(10,*) OUTPUT_MEASUREMENT_FILES
  read(10,*) COMPUTE_ADJOINT_SOURCE
  read(10,*) TSHIFT_MIN, TSHIFT_MAX
  read(10,*) DLNA_MIN, DLNA_MAX
  read(10,*) CC_MIN
  read(10,*) ERROR_TYPE
  read(10,*) DT_SIGMA_MIN
  read(10,*) DLNA_SIGMA_MIN
  read(10,*) ITAPER
  read(10,*) WTR,NPI
  read(10,*) DT_FAC
  read(10,*) ERR_FAC
  read(10,*) DT_MAX_SCALE
  read(10,*) NCYCLE_IN_WINDOW
  read(10,*) USE_PHYSICAL_DISPERSION
  close(10)

  imeas = imeas0
  fstart0 = 1./TLONG ; fend0 = 1./TSHORT
  if(myrank == 0) then 

    ! check the read-in values
    print *, 'INPUTS FROM MEASUREMENT.PAR :'
    print *, '  tt, dtt, nn : ',sngl(tt),sngl(dtt),nn
    print *, '  imeas : ',imeas
    print *, '  chan : ',chan
    print *, '  TLONG, TSHORT : ',sngl(TLONG), sngl(TSHORT)
    
    print *, '  fstart, fend : ', sngl(fstart0), sngl(fend0)
    print *, '  RUN_BANDPASS : ',RUN_BANDPASS
    print *, '  DISPLAY_DETAILS : ',DISPLAY_DETAILS
    print *, '  OUTPUT_MEASUREMENT_FILES : ',OUTPUT_MEASUREMENT_FILES
    print *, '  COMPUTE_ADJOINT_SOURCE : ',COMPUTE_ADJOINT_SOURCE
    print *, '  TSHIFT_MIN, TSHIFT_MAX : ',sngl(TSHIFT_MIN), sngl(TSHIFT_MAX)
    print *, '  DLNA_MIN, DLNA_MAX : ',sngl(DLNA_MIN), sngl(DLNA_MAX)
    print *, '  CC_MIN : ',sngl(CC_MIN)
    print *, '  ERROR_TYPE : ',ERROR_TYPE
    print *, '  DT_SIGMA_MIN : ',sngl(DT_SIGMA_MIN)
    print *, '  DLNA_SIGMA_MIN : ',sngl(DLNA_SIGMA_MIN)
    print *, '  ITAPER : ',ITAPER
    print *, '  WTR, NPI : ',sngl(WTR),NPI
    print *, '  DT_FAC : ',sngl(DT_FAC)
    print *, '  ERR_FAC : ',sngl(ERR_FAC)
    print *, '  DT_MAX_SCALE : ',sngl(DT_MAX_SCALE)
    print *, '  NCYCLE_IN_WINDOW : ',NCYCLE_IN_WINDOW
  endif
  !stop 'checking PAR file input'

! old format way..
!  open(10,file='MEASUREMENT.PAR',status='old',iostat=ios)
!  read(10,'(a)') out_dir
!  read(10,*) is_mtm0
!  read(10,*) wtr,npi
!  read(10,*) iker0
!  read(10,*) RUN_BANDPASS
!  read(10,*) TLONG, TSHORT
!  read(10,*) tt,dtt,nn
!  read(10,*) DISPLAY_DETAILS
!  read(10,*) OUTPUT_MEASUREMENT_FILES
!  read(10,*) INCLUDE_ERROR
!  read(10,*) DT_FAC
!  read(10,*) ERR_FAC
!  read(10,*) DT_MAX_SCALE
!  read(10,*) NCYCLE_IN_WINDOW
!  read(10,*) BEFORE_QUALITY, AFTER_QUALITY
!  read(10,*) BEFORE_TSHIFT, AFTER_TSHIFT
!  read(10,*) DT_SIGMA_MIN, DLNA_SIGMA_MIN
!  close(10)
!
!  out_dir = adjustl(out_dir)
!  iker = iker0
!  is_mtm = is_mtm0
!
!  ! check the read-in values
!  print *, 'INPUTS FROM MEASUREMENT.PAR :'
!  print *, '  is_mtm : ',is_mtm
!  print *, '  wtr, npi : ',wtr,npi
!  print *, '  iker : ',iker
!  print *, '  RUN_BANDPASS :',RUN_BANDPASS
!  print *, '  TLONG, TSHORT : ',TLONG, TSHORT
!  fstart0 = 1./TLONG ; fend0 = 1./TSHORT
!  print *, '  fstart, fend :', fstart0, fend0
!  print *, '  tt, dtt, nn : ',tt,dtt,nn
!  print *, '  out_dir : ',trim(out_dir)
!  print *, '  DISPLAY_DETAILS :',DISPLAY_DETAILS
!  print *, '  OUTPUT_MEASUREMENT_FILES :',OUTPUT_MEASUREMENT_FILES
!  print *, '  INCLUDE_ERROR :',INCLUDE_ERROR
!  print *, '  DT_FAC :',DT_FAC
!  print *, '  ERR_FAC :',ERR_FAC
!  print *, '  DT_MAX_SCALE :',DT_MAX_SCALE
!  print *, '  NCYCLE_IN_WINDOW :',NCYCLE_IN_WINDOW
!  print *, '  BEFORE_QUALITY, AFTER_QUALITY :',BEFORE_QUALITY, AFTER_QUALITY
!  print *, '  BEFORE_TSHIFT, AFTER_TSHIFT :',BEFORE_TSHIFT, AFTER_TSHIFT
!  print *, '  DT_SIGMA_MIN, DLNA_SIGMA_MIN :',DT_SIGMA_MIN, DLNA_SIGMA_MIN
!  !stop 'checking PAR file input'
! apply filter (this should EXACTLY match the filter used in the windowing code)
!trbdndw = 0.3
!a = 30.
!iord = 4
!passes = 2

  ! ray density
  if ( DO_RAY_DENSITY_SOURCE ) ERROR_TYPE = 0

  ! assign additional parameters and stop for certain inconsistencies
  if (fstart0 >= fend0) &
      stop 'Check input frequency range of the signal'

  if (nn > NDIM) &
      stop 'Error: Change interpolation nn or NDIM'

  ! LQY: what happens if imeas = 7/8, and itaper = 2,3, is that permitted?
  ! LQY: how about imeas = 1/2, itaper = 1,2,3 matters?

  if ( imeas == 1 .or. imeas == 2 ) then ! waveforms
    is_mtm0 = 0
  else if ( imeas >= 3 .and. imeas <= 6 ) then ! CC Dt/DlnA
    ! for CC kernels, ITAPER must be a single taper (2 or 3)
    if ( ITAPER == 1 ) stop  'Error: Change ITAPER to 2/3 for CC measurements'
    is_mtm0 = ITAPER     ! 2 or 3 for CC adjoint sources
  else if ( imeas == 7 .or. imeas == 8 ) then
    is_mtm0 = 1          ! multitaper required for MT adjoint source
  else
    stop 'Error: imeas must by 1-8'
  endif

  is_mtm = is_mtm0
  if(myrank == 0) then 
    print *, '  is_mtm :',is_mtm
    print *, ' '
  endif

end subroutine read_par_file_mpi

program measure_adj_mpi

  use measure_adj_mod, only : measure_adj
  use ma_sub,only : NCHI,NDIM,TOL,imeas
  use sacio,only :sacio_readhead,sacio_readsac,sachead
  use mpi

  implicit none

  integer :: npts,nn,npairs,npt1,npt2 
  double precision :: fstart0,fend0,tt,dtt,t0,dt,t01,dt1,t02,dt2  
  CHARACTER(len=10) :: chan,chan_dat,chan_syn
  character(len=512) :: datafile,synfile,line
  character(len=10),dimension(:),allocatable :: netwk,stnm,channel
  character(len=512),dimension(:),allocatable :: dfile,sfile
  double precision,dimension(:,:),allocatable :: window_chi,temp_window
  double precision,dimension(NDIM) :: obs,syn
  double precision,dimension(:),allocatable :: tstart_all,tend_all,&
                                                tr_chi_all,am_chi_all,&
                                                temp_misfit 
  real,dimension(:),allocatable :: temp_data

  type(sachead) :: head

  ! mpi 
  double precision :: tstart,tend
  integer :: nprocs,myrank,ier,ios,i,ipair

  !init mpi
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank,ier)

  ! read par file
  call read_par_file_mpi(myrank,fstart0,fend0,tt,dtt,nn,chan)

  ! input file: MEASUREMENT.WINDOWS
  open(11,file='MEASUREMENT.WINDOWS',status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening input file: MEASUREMENT WINDOWS'

  read(11,*,iostat=ios) npairs
  if (ios /= 0) stop 'Error reading number of pairs of data/syn'
  !print *, 'reading in the data and synthetics...'

  ! allocate space 
  allocate(netwk(npairs),stnm(npairs),window_chi(NCHI,npairs),&
          temp_window(NCHI,npairs),dfile(npairs),sfile(npairs),&
          tstart_all(npairs),tend_all(npairs),&
          tr_chi_all(npairs),am_chi_all(npairs),channel(npairs),&
          temp_misfit(npairs))
  window_chi(:,:) = 0.
  tr_chi_all(:) = 0
  am_chi_all(:) = 0.

  ! loop each pair to read useful information
  do ipair = 1,npairs
    read(11,'(a)',iostat=ios) datafile
    read(11,'(a)',iostat=ios) synfile
    dfile(ipair) = datafile
    sfile(ipair) = synfile

    !print*,trim(datafile)
    !print*,trim(synfile)
    call sacio_readhead(trim(synfile),head,ier)
    netwk(ipair) = trim(head%knetwk)
    stnm(ipair) = trim(head%kstnm)
    channel(ipair) = head%kcmpnm(1:3)

    if(ipair == 1) then ! read useful information
      call sacio_readhead(trim(synfile),head,ier)
      npt1 = head %npts 
      dt1 = head % delta 
      t01 = head%b 

      call sacio_readhead(trim(synfile),head,ier)
      npt2 = head %npts 
      dt2 = head % delta 
      t02 = head%b 
      if (max(npt1,npt2) > NDIM) &
           stop 'Error: Too many number of points in data or syn'
      npts = min(npt1,npt2)

      if (abs(dt1-dt2) > TOL) stop 'Error: check if dt match'
      dt = dt1

      if (abs(t01-t02) > dt)  stop 'Check if t0 match'
      t0 = t01
    endif

    read(11,*)line
    read(11,*) tstart_all(ipair),tend_all(ipair)
  enddo
  close(11)

  ! loop every pair to compute misfits
  do i = myrank,npairs - 1,nprocs 
    ipair = i + 1
    datafile = dfile(ipair)
    synfile = sfile(ipair)
    call sacio_readsac(trim(datafile),head,temp_data,ier)
    obs(1:head%npts) = temp_data(:)
    chan_dat = head %kcmpnm

    call sacio_readsac(trim(synfile),head,temp_data,ier)
    syn(1:head%npts) = temp_data(:)
    tstart = tstart_all(ipair)
    tend = tend_all(ipair)
    chan_syn = head %kcmpnm

    ! compute misfit
    call measure_adj(obs,syn,npts,t0,dt,fstart0,fend0,&
                  tstart,tend,tt,dtt,nn,chan_dat,chan_syn,chan,1,datafile,synfile,&
                  netwk(ipair),stnm(ipair),window_chi(:,ipair),tr_chi_all(ipair),&
                  am_chi_all(ipair))
  enddo

  ! get all window_chi
  temp_window(:,:) = window_chi(:,:)
  call MPI_REDUCE(temp_window,window_chi,NCHI*npairs,MPI_DOUBLE,&
                  MPI_SUM,0,MPI_COMM_WORLD,ier)
  
  temp_misfit(:) = tr_chi_all(:)
  call MPI_REDUCE(temp_misfit,tr_chi_all,npairs,MPI_DOUBLE,&
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  temp_misfit(:) = am_chi_all(:)
  call MPI_REDUCE(temp_misfit,am_chi_all,npairs,MPI_DOUBLE,&
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  if(myrank == 0) then 
    ! ! output files
    ! open(12,file='window_index',status='unknown',iostat=ios)
    open(13,file='window_chi',status='unknown',iostat=ios)
    do ipair = 1,npairs

      synfile = trim(stnm(ipair)) //'.' // trim(netwk(ipair))//'.' // trim(channel(ipair))
      synfile = trim(adjustl(synfile))
      write(13,'(a,x,a,x,a,x,a,x,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
          adjustl(trim(synfile)),trim(stnm(ipair)),trim(netwk(ipair)),trim(channel(ipair)),1,imeas, &
          tstart_all(ipair),tend_all(ipair),window_chi(:,ipair),&
          tr_chi_all(ipair),am_chi_all(ipair),0.,0.

      ! print on screen
      print*,trim(sfile(ipair))
      print*,'Measurement window No. ', 1
      print*,'start and end time of window:',tstart_all(ipair),tend_all(ipair)
      print*,'adjoint source and chi value for imeas =', imeas
      print*,window_chi(7,ipair)
      print *, '     tr_chi = ', sngl(tr_chi_all(ipair)), '  am_chi = ', sngl(am_chi_all(ipair))
      print*,''
    enddo
    close(13)
  endif 

  ! free memory
  deallocate(netwk,stnm,window_chi,temp_window,tstart_all,tend_all,&
            dfile,sfile,tr_chi_all,am_chi_all,channel,temp_misfit)
  if(allocated(temp_data)) deallocate(temp_data)
  call MPI_FINALIZE(ier)


end program measure_adj_mpi
