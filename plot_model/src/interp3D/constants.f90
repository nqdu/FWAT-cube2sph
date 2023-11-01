module constants 

! constants defined
integer,PARAMETER   :: NGLLX = 5,NGLLY = NGLLX, NGLLZ = NGLLX
integer,parameter   :: NGLL3 = NGLLX * NGLLY * NGLLZ  
integer,PARAMETER   :: CUSTOM_REAL = 4,DOUBLE_REAL = 8
integer,PARAMETER   :: UTM_PROJECTION_ZONE = 10
integer,parameter   :: MAX_STR_LEN = 256
integer,PARAMETER   :: IIN = 10,IOUT = 20
double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9
double precision, parameter :: ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0, HALF = 0.5d0
double precision, parameter :: GAUSSALPHA=0.d0,GAUSSBETA=0.d0
integer, parameter :: MIDX = (NGLLX+1)/2
integer, parameter :: MIDY = (NGLLY+1)/2
integer, parameter :: MIDZ = (NGLLZ+1)/2
integer, parameter :: NGNOD_EIGHT_CORNERS = 8
integer, parameter :: NDIM = 3
integer, parameter :: NUM_ITER = 5
integer, parameter :: IDOMAIN_ACOUSTIC    = 1
integer, parameter :: IDOMAIN_ELASTIC     = 2
integer, parameter :: IDOMAIN_POROELASTIC = 3

end module constants