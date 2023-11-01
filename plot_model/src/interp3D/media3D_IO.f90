module read_mesh3D
use interp3D_par
use constants
contains
  subroutine read_external_mesh_for_initialize(MESH_PATH)
    implicit none
    character(len=MAX_STR_LEN)  :: MESH_PATH, filename
    !read nspec and nglob
    write(filename,'(a,a,i6.6,a)') trim(MESH_PATH),'/proc',rank,'_external_mesh.bin'
    open(IIN,file=trim(filename),form='unformatted')
      read(IIN) NSPEC_AB ; read(IIN) NGLOB_AB ; read(IIN) NSPEC_IRREGULAR
    close(IIN)
    allocate(irregular_element_number(NSPEC_AB),stat=ier)
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (NSPEC_IRREGULAR > 0) then
      allocate(xix(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 1'
      allocate(xiy(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 2'
      allocate(xiz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 3'
      allocate(etax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 4'
      allocate(etay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 5'
      allocate(etaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 6'
      allocate(gammax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 7'
      allocate(gammay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 8'
      allocate(gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 9'
      allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) stop 'cannot allocate 10'
    else
      allocate(xix(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 1'
      allocate(xiy(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 2'
      allocate(xiz(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 3'
      allocate(etax(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 4'
      allocate(etay(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 5'
      allocate(etaz(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 6'
      allocate(gammax(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 7'
      allocate(gammay(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 7'
      allocate(gammaz(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 9'
      allocate(jacobian(1,1,1,1),stat=ier)
      if (ier /= 0) stop 'cannot allocate 10'
    endif

    ! mesh node locations
    allocate(xstore(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 11'
    allocate(ystore(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 12'
    allocate(zstore(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 13'

    ! material properties
    allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 14'
    allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 15'
    ! material flags
    allocate(ispec_is_acoustic(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 16'
    allocate(ispec_is_elastic(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 17'
    allocate(ispec_is_poroelastic(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'cannot allocate 18'
    ispec_is_acoustic(:) = .false.
    ispec_is_elastic(:) = .false.
    ispec_is_poroelastic(:) = .false.
  end subroutine read_external_mesh_for_initialize
  
  subroutine read_external_mesh(MESH_PATH)
    implicit none
    character(len=MAX_STR_LEN)  :: MESH_PATH, filename
    !read nspec and nglob
    write(filename,'(a,a,i6.6,a)') trim(MESH_PATH),'/proc',rank,'_external_mesh.bin'
    open(IIN,file=trim(filename),form='unformatted')
      read(IIN) NSPEC_AB
      read(IIN) NGLOB_AB
      read(IIN) NSPEC_IRREGULAR

      read(IIN) ibool

      read(IIN) xstore
      read(IIN) ystore
      read(IIN) zstore

      read(IIN) irregular_element_number
      read(IIN) xix_regular
      read(IIN) jacobian_regular

      read(IIN) xix
      read(IIN) xiy
      read(IIN) xiz
      read(IIN) etax
      read(IIN) etay
      read(IIN) etaz
      read(IIN) gammax
      read(IIN) gammay
      read(IIN) gammaz
      read(IIN) jacobian

      read(IIN) kappastore
      read(IIN) mustore

      read(IIN) ispec_is_acoustic
      read(IIN) ispec_is_elastic
      read(IIN) ispec_is_poroelastic
    close(IIN)
    !print *,'xix',xix(1,1,1,1),xix_regular,'rank=',rank
  end subroutine read_external_mesh

end module read_mesh3D
    
module media_IO
  
  contains
  subroutine read_coor_get_lines(filename)
    use constants
    use interp3D_par
    implicit none
  
    character(len=MAX_STR_LEN) :: filename 
    ! local
    character(len=MAX_STR_LEN) :: line 

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

  subroutine read_coordinate(COOR_PATH,SUPRESS_UTM)
    use constants
    use interp3D_par
    implicit none
    character(len=MAX_STR_LEN) :: COOR_PATH,SUPRESS_UTM
    real(kind=8)  :: rlon,rlat,rx,ry
    integer :: ipts
    allocate(model_in(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    allocate(xin(npts),stat=ier)
    allocate(yin(npts),stat=ier)
    allocate(zin(npts),stat=ier)
    allocate(model_out(npts),stat=ier)
    allocate(x_xi(npts),stat=ier)
    allocate(y_eta(npts),stat=ier)
    allocate(z_gamma(npts),stat=ier)
    allocate(dist_xyz_newxyz(npts),stat=ier)
    open(IIN,file=COOR_PATH)
      do ipts=1,npts
        read(IIN,*)rlon,rlat,zin(ipts)
        if(SUPRESS_UTM=='false')then
          call utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,0) ! geo2utm
          xin(ipts) = real(rx,kind=CUSTOM_REAL)
          yin(ipts) = real(ry,kind=CUSTOM_REAL)
        else
          xin(ipts) = real(rlon,kind=CUSTOM_REAL)
          yin(ipts) = real(rlat,kind=CUSTOM_REAL)
        endif
      enddo
    close(IIN)

  end subroutine read_coordinate
  
  subroutine read_binORdat_media(binORdat,MODEL_PATH,PAR_NAME)
    use constants
    use interp3D_par
    
    implicit none
    character(len=MAX_STR_LEN)  :: MODEL_PATH, filename, binORdat, PAR_NAME
    integer :: ispec,i,j,k
    if(binORdat=='bin')then
      write(filename,'(a,a,i6.6,a,a,a)') trim(MODEL_PATH),'/proc',rank,'_',trim(PAR_NAME),'.bin'
      open(IIN,file=trim(filename),form='unformatted',iostat=ier)
      if(ier /=0 ) stop 'error readind vpfile '; read(IIN) model_in; close(IIN)
      if(rank == 0)  then 
        print*,'finished reading parameters '//trim(filename)
      endif
    elseif(binORdat=='dat')then
      write(filename,'(a,i5.5,a,i5.5,a)')trim(MODEL_PATH)//'/'//trim(PAR_NAME)//'_',rank,'.dat'
      open(IIN,file=trim(filename),iostat=ier,form='formatted')
      if(ier /=0 ) stop 'error readind vpfile ' 
      do ispec = 1,NSPEC_AB
       do k=1,NGLLZ  
        do j=1,NGLLY 
         do i = 1,NGLLX 
          read(IIN,*) model_in(i,j,k,ispec)
         enddo
        enddo
       enddo
      enddo
      close(IIN)
      if(rank == 0)  then 
        print*,'finished reading parameters '//trim(filename)
      endif
    else
      if(rank == 0)  then 
        stop 'wrong parameter for datatype, either bin for dat '
      endif
    endif 
    end subroutine read_binORdat_media

    subroutine write_interped_media(OUT_NAME,SUPRESS_UTM)
    use constants
    use interp3D_par
    use ieee_arithmetic,only : ieee_value,ieee_quiet_nan
    implicit none
    character(len=MAX_STR_LEN) :: OUT_NAME,SUPRESS_UTM
    real(kind=8)  :: rlon,rlat,rx,ry
    integer :: ipts
    where(model_out == 0.0)
      model_out = ieee_value(model_out,ieee_quiet_nan)
    endwhere
    open(IIN,file=OUT_NAME)
      do ipts=1,npts
        rx=xin(ipts)
        ry=yin(ipts)
        if(SUPRESS_UTM=='false')then
          call utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,1) ! geo2utm
          rlon = real(rlon,kind=CUSTOM_REAL)
          rlat = real(rlat,kind=CUSTOM_REAL)
          write(IIN,*)rlon,rlat,zin(ipts) * 0.001,model_out(ipts)
        else
          write(IIN,*)xin(ipts),yin(ipts),zin(ipts) * 0.001,model_out(ipts)
        endif
      enddo
    close(IIN)
    end subroutine write_interped_media


end module media_IO