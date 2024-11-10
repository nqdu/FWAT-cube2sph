subroutine seis_pca(stf_collect,nrec,npts,stf) bind(c,name='seis_pca')
  use iso_c_binding
  use spanlib,only : sl_pca
  implicit none
  ! variables
  integer(c_int),value,intent(in) :: nrec,npts 
  real(c_float),intent(in) :: stf_collect(npts,nrec)
  real(c_float),intent(inout) :: stf(npts)

  ! local
  real(c_float),allocatable :: xeof(:,:),pc(:,:),pcT(:,:),pdm(:,:),ev(:),ff(:,:)
  integer :: nkeep

  ! allocate space 
  nkeep = nrec
  allocate(xeof(nkeep,nkeep),ff(nrec,npts))
  allocate(pc(npts,nkeep))
  allocate(pcT(nkeep,npts))
  allocate(pdm(npts,nkeep))
  allocate(ev(nkeep))
  ff = transpose(stf_collect)
  call sl_pca(ff, nkeep, xeof, pc, ev) 

  ! copy results to stf
  stf(:) = pc(:,1)

  ! deallocate
  deallocate(xeof,pc,pcT,pdm,ev,ff)
end subroutine 