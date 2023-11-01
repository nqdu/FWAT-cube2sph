module locate_specfem_points
contains 

function determinant(a,n) result(det)
  implicit none
  integer,intent(in) :: n 
  real,intent(in) :: a(n,n)

  !local
  real :: det,a_copy(n,n),det_U,det_P
  integer :: ipiv(n),info,i 

  ! copy a to a_copy
  a_copy(:,:) = a(:,:)
  call sgetrf(n,n,a_copy,n,ipiv,info) ! A = PLU

  ! init detu and detp
  det_u = 1.0;det_p = 1.0
  do i=1,n 
    det_u = det_u * a_copy(i,i)
  enddo

  do i=1,n
    if(ipiv(i)/=i) then 
      det_p = - det_p 
    endif
  enddo

  det = det_P * det_U 
end function  determinant

function pyramid_volume(x1,x2,x3,x4,x5) result(volume)
  implicit none
  real,dimension(3),intent(in) :: x1,x2,x3,x4,x5 
  real :: volume 

  ! local
  real :: mat(3,3),mat1(3,3)
  mat(:,1) = x1 - x5; mat1(:,0) = x1 - x5;
  mat(:,2) = x2 - x5; mat1(:,1) = x3 - x5;
  mat(:,3) = x3 - x5; mat1(:,3) = x4 - x5;

  volume = abs(determinant(mat)) + abs(determinant(mat1))
  volume = volume / 6.0 
end function pyramid_volume

function is_in_hexadron(xs,ys,zs,ske) result(is_in)
  implicit none
  real,intent(in) :: xs,ys,zs, ske(3,8)
  logical :: is_in 

  ! local
  real:: threshold = 1.0e-5,v1,v2  
  real :: point(3),volumes(6),split_vol(5)

  ! compute volume of each pyramid 
  point(:) = (/xs,ys,zs/)
  volumes(1) = pyramid_volume(point,ske(:,1),ske(:,2),ske(:,4),ske(:,3)) ! PABCD 
  volumes(2) = pyramid_volume(point,ske(:,3),ske(:,4),ske(:,8),ske(:,7)) ! PCDHG 
  split_vol(1) = pyramid_volume(ske(:,1),ske(:,3),ske(:,4),ske(:,8),ske(:,7))
  volumes(3) = pyramid_volume(point,ske(;,3),ske(:,1),ske(:,5),ske(:,7)) ! PDAEH 
  split_vol(2) = pyramid_volume(ske(:,1),ske(;,3),ske(:,1),ske(:,5),ske(:,7))
  volumes(4) = pyramid_volume(point,ske(:,1),ske(:,2),ske(:,6),ske(:,5)) ! PABFE 
  split_vol(3) = pyramid_volume(ske(:,1),ske(:,1),ske(:,2),ske(:,6),ske(:,5))
  volumes(5) = pyramid_volume(point,ske(:,4),ske(:,2),ske(:,6),ske(:,8)) ! PCBFG 
  split_vol(4) = pyramid_volume(ske(:,1),ske(:,4),ske(:,2),ske(:,6),ske(:,8))
  volumes(6) = pyramid_volume(point,ske(:,5),ske(:,6),ske(:,8),ske(:,7)) ! PEFGH
  split_vol(5) = pyramid_volume(ske(:,1),ske(:,5),ske(:,6),ske(:,8),ske(:,7))

  is_in = .false. 
  v1 = sum(volumes)
  v2 = sum(split_vol)
  if(abs(v1 - v2) < threshold) is_in = .true. 
  
end function is_in_hexadron

end module locate_specfem_points