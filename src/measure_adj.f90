subroutine measure_adj_fwat(obs_in,syn_in,npts,t0,dt,&
                          verbose,tstart,tend,tt,dtt,nn,&
                          chan_c, imeas0_in,TLONG_in, TSHORT_in, &
                          RUN_BANDPASS_in, DISPLAY_DETAILS_in, OUTPUT_MEASUREMENT_FILES_in, &
                          COMPUTE_ADJOINT_SOURCE_in, TSHIFT_MIN_in, TSHIFT_MAX_in, &
                          DLNA_MIN_in, DLNA_MAX_in, CC_MIN_in, ERROR_TYPE_in, &
                          DT_SIGMA_MIN_in, DLNA_SIGMA_MIN_in, ITAPER_in, &
                          WTR_in, NPI_in, DT_FAC_in, ERR_FAC_in, DT_MAX_SCALE_in, &
                          NCYCLE_IN_WINDOW_in, USE_PHYSICAL_DISPERSION_in,&
                          window_chi,tr_chi,am_chi,adj_out) bind(C,name='measure_adj_fwat_')
  use iso_c_binding
  use ma_constants
  use ma_variables
  use ma_sub2
  use ascii_rw
  use measure_adj_mod,only: measure_adj
  implicit none
  

  integer(c_int),value,intent(in) :: npts,nn
  real(c_double),value,intent(in) :: tstart,tend,t0,dt,tt,dtt
  real(c_double),value, intent(in) :: TLONG_in, TSHORT_in
  real(c_double),value, intent(in) :: TSHIFT_MIN_in, TSHIFT_MAX_in
  real(c_double),value, intent(in) :: DLNA_MIN_in, DLNA_MAX_in
  real(c_double),value, intent(in) :: CC_MIN_in, DT_SIGMA_MIN_in, DLNA_SIGMA_MIN_in
  real(c_double),value, intent(in) :: WTR_in, NPI_in
  real(c_double),value, intent(in) :: DT_FAC_in, ERR_FAC_in, DT_MAX_SCALE_in
  real(c_double),value, intent(in) :: NCYCLE_IN_WINDOW_in
  integer(c_int),value, intent(in) :: imeas0_in, ERROR_TYPE_in, ITAPER_in
  logical(c_bool),value, intent(in) :: RUN_BANDPASS_in, DISPLAY_DETAILS_in, OUTPUT_MEASUREMENT_FILES_in
  logical(c_bool),value, intent(in) :: COMPUTE_ADJOINT_SOURCE_in, USE_PHYSICAL_DISPERSION_in,verbose 
  character(kind=c_char) :: chan_c(*)
  real(c_double),intent(in) :: obs_in(npts),syn_in(npts)
  real(c_double),intent(inout) :: window_chi(NCHI),tr_chi,am_chi,adj_out(nn)

  ! local
  double precision :: obs(NDIM),syn(NDIM),fstart0,fend0 
  integer :: j,len_chan_c
  CHARACTER(len=10) :: chan,chan_dat,chan_syn,net,stnm
  character(len=128) :: datafile,synfile

  ! copy dat 
  obs(1:npts) = obs_in(:)
  syn(1:npts) = syn_in(:)

  ! copy to global vars
  ! Assignments
  imeas0 = imeas0_in
  TLONG = TLONG_in
  TSHORT = TSHORT_in
  RUN_BANDPASS = RUN_BANDPASS_in
  DISPLAY_DETAILS = DISPLAY_DETAILS_in
  OUTPUT_MEASUREMENT_FILES = OUTPUT_MEASUREMENT_FILES_in
  COMPUTE_ADJOINT_SOURCE = COMPUTE_ADJOINT_SOURCE_in
  USE_PHYSICAL_DISPERSION = USE_PHYSICAL_DISPERSION_in
  TSHIFT_MIN = TSHIFT_MIN_in
  TSHIFT_MAX = TSHIFT_MAX_in
  DLNA_MIN = DLNA_MIN_in
  DLNA_MAX = DLNA_MAX_in
  CC_MIN = CC_MIN_in
  ERROR_TYPE = ERROR_TYPE_in
  DT_SIGMA_MIN = DT_SIGMA_MIN_in
  DLNA_SIGMA_MIN = DLNA_SIGMA_MIN_in
  ITAPER = ITAPER_in
  WTR = WTR_in
  NPI = NPI_in
  DT_FAC = DT_FAC_in
  ERR_FAC = ERR_FAC_in
  DT_MAX_SCALE = DT_MAX_SCALE_in
  NCYCLE_IN_WINDOW = NCYCLE_IN_WINDOW_in
  OUT_DIR = 'OUTPUT_FILES'   ! default

  ! get chan 
  do j = 1,10 
    if (chan_c(j) == c_null_char) then
      len_chan_c = j - 1
      exit
    end if
  end do
  chan = transfer(chan_c(1:len_chan_c), chan)

  imeas = imeas0
  fstart0 = 1./TLONG ; fend0 = 1./TSHORT
  if(verbose) then 

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
  if(verbose) then 
    print *, '  is_mtm :',is_mtm
    print *, ' '
  endif
  
  chan_dat = chan 
  chan_syn = chan 
  datafile = chan 
  synfile = chan
  net = chan 
  stnm = chan
  call measure_adj(obs,syn,npts,t0,dt,fstart0,fend0,tstart,tend,&
                    tt,dtt,nn,chan_dat,chan_syn,chan,1,datafile,synfile,&
                    net,stnm,window_chi,tr_chi,am_chi,adj_out)

end subroutine measure_adj_fwat
