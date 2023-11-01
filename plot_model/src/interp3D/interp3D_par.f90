module interp3D_par
  use constants
! par for mpi
  integer ier, rank, nproc
! copy from specfem3D_par

 integer :: NGNOD=8
! number of spectral element and global points
  integer :: NSPEC_AB, NGLOB_AB

! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore

  integer :: NSPEC_IRREGULAR
  integer, dimension(:), allocatable :: irregular_element_number
  real(kind=CUSTOM_REAL) :: xix_regular,jacobian_regular

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

! material properties
  ! isotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore,mustore

! density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore
! materials
  logical, dimension(:), allocatable :: ispec_is_poroelastic
  logical, dimension(:), allocatable :: ispec_is_elastic
  logical, dimension(:), allocatable :: ispec_is_acoustic

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprime_yyT,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprime_zzT,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

! Lagrange interpolators at receiver(for this code, they are at points for interp)
  double precision, dimension(:), allocatable :: hxi,heta,hpxi,hpeta,hgamma,hpgamma
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: hxi_store,heta_store,hgamma_store, points_location_store
  integer, dimension(:), allocatable:: ispec_selected_store
! parameters for locate a points in mesh
! number of points for interp3D
  integer :: npts 
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: model_in
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xin,yin,zin,xin_found,yin_found,zin_found,model_out
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: x_xi,y_eta,z_gamma,dist_xyz_newxyz

end module interp3D_par
