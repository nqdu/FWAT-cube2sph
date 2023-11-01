module locate_points

! use interp3D_par
use constants,only : CUSTOM_REAL
implicit none
real(kind=CUSTOM_REAL) :: elemsize_max_glob,elemsize_min_glob

contains

subroutine setup_GLL_points()
  use interp3D_par
  use constants
  implicit none
  integer :: j,i,k
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                                  hprime_xx,hprime_yy,hprime_zz, &
                                  hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                                  wgllwgll_xy,wgllwgll_xz,wgllwgll_yz)
  ! define transpose of derivation matrix
  do j = 1,NGLLY
    do i = 1,NGLLX
      hprime_xxT(j,i) = hprime_xx(i,j)
      hprime_yyT(j,i) = hprime_yy(i,j)
      hprime_zzT(j,i) = hprime_zz(i,j)
      hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
    enddo
  enddo

  ! define a 3D extension in order to be able to force vectorization in the compute_forces routines
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        wgllwgll_yz_3D(i,j,k) = wgllwgll_yz(j,k)
        wgllwgll_xz_3D(i,j,k) = wgllwgll_xz(i,k)
        wgllwgll_xy_3D(i,j,k) = wgllwgll_xy(i,j)
      enddo
    enddo
  enddo

  ! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxi(NGLLX),stat=ier)
  allocate(hpxi(NGLLX),stat=ier)
  allocate(heta(NGLLY),stat=ier)
  allocate(hpeta(NGLLY),stat=ier)
  allocate(hgamma(NGLLZ),stat=ier)
  allocate(hpgamma(NGLLZ),stat=ier)
  
end subroutine setup_GLL_points

subroutine get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                        NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

! calculates the min/max size of an edge of the specified element (ispec);
! we purposely do not include the distance along the diagonals of the element, only the size of its edges.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL

  implicit none

  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: dist
  real(kind=CUSTOM_REAL) :: x1,x2,y1,y2,z1,z2
  integer :: i,j,k
  integer :: i1,i2,j1,j2,k1,k2
  integer :: iglob1,iglob2

  ! initializes
  elemsize_min = HUGEVAL
  elemsize_max = -HUGEVAL

  ! loops over the four edges that are along X
  i1 = 1
  i2 = NGLLX
  do k = 1, NGLLZ, NGLLZ-1
    do j = 1, NGLLY, NGLLY-1
      iglob1 = ibool(i1,j,k,ispec)
      iglob2 = ibool(i2,j,k,ispec)

      x1 = xstore(iglob1)
      y1 = ystore(iglob1)
      z1 = zstore(iglob1)

      x2 = xstore(iglob2)
      y2 = ystore(iglob2)
      z2 = zstore(iglob2)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

      if (dist < elemsize_min) elemsize_min = dist
      if (dist > elemsize_max) elemsize_max = dist

    enddo
  enddo

  ! loops over the four edges that are along Y
  j1 = 1
  j2 = NGLLY
  do k = 1, NGLLZ, NGLLZ-1
    do i = 1, NGLLX, NGLLX-1
      iglob1 = ibool(i,j1,k,ispec)
      iglob2 = ibool(i,j2,k,ispec)

      x1 = xstore(iglob1)
      y1 = ystore(iglob1)
      z1 = zstore(iglob1)

      x2 = xstore(iglob2)
      y2 = ystore(iglob2)
      z2 = zstore(iglob2)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

      if (dist < elemsize_min) elemsize_min = dist
      if (dist > elemsize_max) elemsize_max = dist
    enddo
  enddo

  ! loops over the four edges that are along Z
  k1 = 1
  k2 = NGLLZ
  do j = 1, NGLLY, NGLLY-1
    do i = 1, NGLLX, NGLLX-1
      iglob1 = ibool(i,j,k1,ispec)
      iglob2 = ibool(i,j,k2,ispec)

      x1 = xstore(iglob1)
      y1 = ystore(iglob1)
      z1 = zstore(iglob1)

      x2 = xstore(iglob2)
      y2 = ystore(iglob2)
      z2 = zstore(iglob2)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

      if (dist < elemsize_min) elemsize_min = dist
      if (dist > elemsize_max) elemsize_max = dist
    enddo
  enddo
  elemsize_min = sqrt( elemsize_min )
  elemsize_max = sqrt( elemsize_max )

end subroutine get_elem_minmaxsize

subroutine get_elem_minmaxsize_glob()
  use constants
  use interp3D_par
  use mpi_f08
  implicit none

  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max
  integer :: ispec
  ! initializes
  elemsize_min_glob = HUGEVAL
  elemsize_max_glob = -HUGEVAL
  do ispec=1,NSPEC_AB
    call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                        NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)
    elemsize_min_glob = min(elemsize_min_glob,elemsize_min)
    elemsize_max_glob = max(elemsize_max_glob,elemsize_max)
  enddo
  elemsize_min = elemsize_min_glob
  elemsize_max = elemsize_max_glob
  
  call MPI_REDUCE(elemsize_min,elemsize_min_glob,1,MPI_REAL,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call MPI_REDUCE(elemsize_max,elemsize_max_glob,1,MPI_REAL,MPI_MIN,0,MPI_COMM_WORLD,ier)
  if(rank==0)print*,'glob elemsize min and max is:',elemsize_min_glob,elemsize_max_glob
end subroutine get_elem_minmaxsize_glob

subroutine locate_points_in_mesh(x_target, y_target, z_target, elemsize_max_glob, & 
                                  ispec_selected, xi_found, eta_found, gamma_found, &
                                  x_found, y_found, z_found, domain, final_distance_squared)
  use constants
  use interp3D_par
  implicit none
  
  double precision,                      intent(in)     :: x_target, y_target, z_target
  real(kind=CUSTOM_REAL),                intent(in)     :: elemsize_max_glob
  double precision,                      intent(out)    :: x_found,  y_found,  z_found
  double precision,                      intent(out)    :: xi_found, eta_found, gamma_found
  integer,                               intent(out)    :: ispec_selected, domain
  double precision,                      intent(out)    :: final_distance_squared

  ! locals
  integer                                               :: iter_loop , ispec, iglob, i, j, k, ia, iax, iay, iaz
  ! location search
  double precision                                      :: maximal_elem_size_squared, dist_squared
  double precision                                      :: distmin_squared
  double precision                                      :: x,y,z
  double precision                                      :: xi,eta,gamma,dx,dy,dz,dxi,deta
  double precision                                      :: xixs,xiys,xizs
  double precision                                      :: etaxs,etays,etazs
  double precision                                      :: gammaxs,gammays,gammazs, dgamma
  integer                                               :: ix_initial_guess, iy_initial_guess, iz_initial_guess
  integer, dimension(NGNOD)                             :: iaddx,iaddy,iaddz
  double precision, dimension(NGNOD)                    :: xelm,yelm,zelm
  double precision :: final_distance_squared_this_element
  logical, dimension(NGLOB_AB) :: flag_topological
  logical :: use_adjacent_elements_search
  integer :: number_of_mesh_elements_for_the_initial_guess
  integer, dimension(50) :: array_of_all_elements_of_ispec_selected   ! 15% faster
  logical, parameter :: USE_SINGLE_PASS = .true.

  maximal_elem_size_squared = (10. * elemsize_max_glob)**2
  ispec_selected   = 1
  ix_initial_guess = 1
  iy_initial_guess = 1
  iz_initial_guess = 1
  distmin_squared = HUGEVAL

  do ispec = 1, NSPEC_AB
    iglob = ibool(MIDX,MIDY,MIDZ,ispec)
    dist_squared = (x_target - dble(xstore(iglob)))**2 &
                 + (y_target - dble(ystore(iglob)))**2 &
                 + (z_target - dble(zstore(iglob)))**2
    
    if (dist_squared > maximal_elem_size_squared) cycle ! exclude elements that are too far from target
    ! find closest GLL point form target
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distance to GLL point
          dist_squared = (x_target - dble(xstore(iglob)))**2 &
                       + (y_target - dble(ystore(iglob)))**2 &
                       + (z_target - dble(zstore(iglob)))**2

          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ispec_selected  = ispec
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
            x_found = xstore(iglob)
            y_found = ystore(iglob)
            z_found = zstore(iglob)
          endif
        enddo
      enddo
    enddo
    
  enddo
  ! checks if initial guess lies on an element boundary
  if ((ix_initial_guess == 1 .or. ix_initial_guess == NGLLX) .or. &
      (iy_initial_guess == 1 .or. iy_initial_guess == NGLLY) .or. &
      (iz_initial_guess == 1 .or. iz_initial_guess == NGLLZ)) then
    ! best guess close to a GLL point on an element boundary
    ! negect if one the surface BinHe
      use_adjacent_elements_search = .true.
  else
    ! best guess lies within the selected element, no need for checking adjacent elements
    use_adjacent_elements_search = .false.
  endif

  ! debug
  !print *,'locate ',use_adjacent_elements_search,ix_initial_guess,iy_initial_guess,iz_initial_guess
  ! flagging corners
  flag_topological(:) = .false.
  ! mark the eight corners of the initial guess element
  flag_topological(ibool(1,1,1,ispec_selected)) = .true.
  flag_topological(ibool(NGLLX,1,1,ispec_selected)) = .true.
  flag_topological(ibool(NGLLX,NGLLY,1,ispec_selected)) = .true.
  flag_topological(ibool(1,NGLLY,1,ispec_selected)) = .true.

  flag_topological(ibool(1,1,NGLLZ,ispec_selected)) = .true.
  flag_topological(ibool(NGLLX,1,NGLLZ,ispec_selected)) = .true.
  flag_topological(ibool(NGLLX,NGLLY,NGLLZ,ispec_selected)) = .true.
  flag_topological(ibool(1,NGLLY,NGLLZ,ispec_selected)) = .true.
  if (use_adjacent_elements_search) then
  if (USE_SINGLE_PASS) then
    ! assumes a maximum of a 100 adjacent elements
    !
    ! first store the initial guess itself
    number_of_mesh_elements_for_the_initial_guess = 1
    array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec_selected

    ! then store all the others
    do ispec = 1, NSPEC_AB

      ! omit selected element since already included
      if (ispec == ispec_selected) cycle

      ! loop on the eight corners only, no need to loop on the rest since we just want to detect adjacency
      do k = 1,NGLLZ,NGLLZ-1
        do j = 1,NGLLY,NGLLY-1
          do i = 1,NGLLX,NGLLX-1
            if (flag_topological(ibool(i,j,k,ispec))) then
              ! this element is in contact with the initial guess
              number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
              ! check
              if (number_of_mesh_elements_for_the_initial_guess > 100) stop 'Error must increase array size in locate_point.f90'

              array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec
              ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
              goto 707
            endif
          enddo
        enddo
      enddo
      707 continue
    enddo

  else
    ! note: this involves 2-passes over all elements and dynamic memory allocation
    !       which can be very slow for large meshes

    ! loop on all the elements to count how many are shared with the initial guess
    number_of_mesh_elements_for_the_initial_guess = 1
    do ispec = 1, NSPEC_AB
      if (ispec == ispec_selected) cycle

      ! loop on the eight corners only, no need to loop on the rest since we just want to detect adjacency
      do k = 1,NGLLZ,NGLLZ-1
        do j = 1,NGLLY,NGLLY-1
          do i = 1,NGLLX,NGLLX-1
            if (flag_topological(ibool(i,j,k,ispec))) then
              ! this element is in contact with the initial guess
              number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
              ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
              goto 700
            endif
          enddo
        enddo
      enddo
      700 continue
    enddo
    ! debug
    !print *,'number_of_mesh_elements_for_the_initial_guess = ',number_of_mesh_elements_for_the_initial_guess

    ! now that we know the number of elements, we can allocate the list of elements and create it
    !allocate(array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess))

    ! first store the initial guess itself
    number_of_mesh_elements_for_the_initial_guess = 1
    array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec_selected

    ! then store all the others
    do ispec = 1, NSPEC_AB

      if (ispec == ispec_selected) cycle

      ! loop on the eight corners only, no need to loop on the rest since we just want to detect adjacency
      do k = 1,NGLLZ,NGLLZ-1
        do j = 1,NGLLY,NGLLY-1
          do i = 1,NGLLX,NGLLX-1
            if (flag_topological(ibool(i,j,k,ispec))) then
              ! this element is in contact with the initial guess
              number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
              array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec
              ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
              goto 800
            endif
          enddo
        enddo
      enddo
      800 continue
    enddo
  endif !end of USE_SINGLE_PASS
  ! debug
  !print *,'adjacent elements ',number_of_mesh_elements_for_the_initial_guess
  else
  ! frees search results
  number_of_mesh_elements_for_the_initial_guess = 1
  array_of_all_elements_of_ispec_selected(1) = ispec_selected
  endif !end of use_adjacent_elements_search

  ! define topology of the control element
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

!! DK DK dec 2017
  final_distance_squared = HUGEVAL

  do i = 1,number_of_mesh_elements_for_the_initial_guess

!! DK DK dec 2017 set initial guess in the middle of the element, since we computed the true one only for the true initial guess
!! DK DK dec 2017 the nonlinear process below will converge anyway
  if (i > 1) then
    ix_initial_guess = (NGLLX+1) / 2
    iy_initial_guess = (NGLLY+1) / 2
    iz_initial_guess = (NGLLZ+1) / 2
  endif

  ispec = array_of_all_elements_of_ispec_selected(i)

  ! general coordinate of initial guess
  xi    = xigll(ix_initial_guess)
  eta   = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  ! define coordinates of the control points of the element
  do ia = 1,NGNOD
    iax = 0
    iay = 0
    iaz = 0
    if (iaddx(ia) == 0) then
        iax = 1
     else if (iaddx(ia) == 1) then
        iax = (NGLLX+1)/2
     else if (iaddx(ia) == 2) then
        iax = NGLLX
     else
        stop 'incorrect value of iaddx'
    endif

    if (iaddy(ia) == 0) then
        iay = 1
    else if (iaddy(ia) == 1) then
      iay = (NGLLY+1)/2
    else if (iaddy(ia) == 2) then
      iay = NGLLY
    else
      stop 'incorrect value of iaddy'
    endif

    if (iaddz(ia) == 0) then
        iaz = 1
    else if (iaddz(ia) == 1) then
      iaz = (NGLLZ+1)/2
    else if (iaddz(ia) == 2) then
        iaz = NGLLZ
    else
        stop 'incorrect value of iaddz'
    endif

    iglob = ibool(iax,iay,iaz,ispec)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))
  enddo

    ! iterate to solve the non linear system
    do iter_loop = 1, NUM_ITER

      ! recompute jacobian for the new point
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                              xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

      ! compute distance to target location
      dx = - (x - x_target)
      dy = - (y - y_target)
      dz = - (z - z_target)

      ! compute increments
      dxi  = xixs*dx + xiys*dy + xizs*dz
      deta = etaxs*dx + etays*dy + etazs*dz
      dgamma = gammaxs*dx + gammays*dy + gammazs*dz

      ! update values
      xi = xi + dxi
      eta = eta + deta
      gamma = gamma + dgamma

      ! impose that we stay in that element
      ! (useful if user gives a point outside the mesh for instance)
      if (xi > 1.01d0) xi     =  1.01d0
      if (xi < -1.01d0) xi     = -1.01d0
      if (eta > 1.01d0) eta    =  1.01d0
      if (eta < -1.01d0) eta    = -1.01d0
      if (gamma > 1.01d0) gamma  =  1.01d0
      if (gamma < -1.01d0) gamma  = -1.01d0

    enddo

    ! compute final coordinates of point found
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

    ! compute final distance squared between asked and found
    final_distance_squared_this_element = (x_target-x)**2 + (y_target-y)**2 + (z_target-z)**2

    ! if we have found an element that gives a shorter distance
    if (final_distance_squared_this_element < final_distance_squared) then
      ! store information about the point found
      ! note: xi/eta/gamma will be in range [-1,1]
      ispec_selected = ispec

      xi_found = xi
      eta_found = eta
      gamma_found = gamma

      x_found = x
      y_found = y
      z_found = z

      !   store final distance squared between asked and found
      final_distance_squared = final_distance_squared_this_element
    endif

! !! DK DK dec 2017
enddo

  ! sets whether acoustic (1) or elastic (2)
  if (ispec_is_acoustic( ispec_selected )) then
    domain = IDOMAIN_ACOUSTIC
  else if (ispec_is_elastic( ispec_selected )) then
    domain = IDOMAIN_ELASTIC
  else if (ispec_is_poroelastic( ispec_selected )) then
    domain = IDOMAIN_POROELASTIC
  else
    domain = 0
  endif

end subroutine locate_points_in_mesh



end module locate_points
