subroutine mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

    ! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3comp_singleB
#else
! cray
!DIR$ INLINEALWAYS mxm5_3comp_singleB
#endif
    
    ! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays
    
      implicit none
      integer,parameter :: CUSTOM_REAL = 4 
    
      integer,intent(in) :: n1,n3
      real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A1,A2,A3
      real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
      real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3
    
      ! local parameters
      integer :: i,j
    
      ! matrix-matrix multiplication
      do j = 1,n3
!DIR$ IVDEP
#if defined __INTEL_COMPILER
!DIR$ SIMD
#endif
        do i = 1,n1
          ! C1(i,j) =  A1(i,1) * B(1,j) &
          !          + A1(i,2) * B(2,j) &
          !          + A1(i,3) * B(3,j) &
          !          + A1(i,4) * B(4,j) &
          !          + A1(i,5) * B(5,j)
    
          ! C2(i,j) =  A2(i,1) * B(1,j) &
          !          + A2(i,2) * B(2,j) &
          !          + A2(i,3) * B(3,j) &
          !          + A2(i,4) * B(4,j) &
          !          + A2(i,5) * B(5,j)
    
          ! C3(i,j) =  A3(i,1) * B(1,j) &
          !          + A3(i,2) * B(2,j) &
          !          + A3(i,3) * B(3,j) &
          !          + A3(i,4) * B(4,j) &
          !          + A3(i,5) * B(5,j)
        C1(i,j) = sum(A1(i,:)*B(:,j))
        C2(i,j) = sum(A2(i,:)*B(:,j))
        C3(i,j) = sum(A3(i,:)*B(:,j))
        enddo
      enddo
    
      end subroutine mxm5_3comp_singleB
    
subroutine mxm5_3comp(A1,A2,A3,B,C1,C2,C3)
    implicit none
    integer,parameter :: CUSTOM_REAL = 4,n1 = 25,n3 = 5
    real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A1,A2,A3
    real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
    real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

      ! local parameters
    integer :: i,j
    
    ! matrix-matrix multiplication
    do j = 1,n3
      do i = 1,n1
        ! C1(i,j) =  A1(i,1) * B(1,j) &
        !          + A1(i,2) * B(2,j) &
        !          + A1(i,3) * B(3,j) &
        !          + A1(i,4) * B(4,j) &
        !          + A1(i,5) * B(5,j)
  
        ! C2(i,j) =  A2(i,1) * B(1,j) &
        !          + A2(i,2) * B(2,j) &
        !          + A2(i,3) * B(3,j) &
        !          + A2(i,4) * B(4,j) &
        !          + A2(i,5) * B(5,j)
  
        ! C3(i,j) =  A3(i,1) * B(1,j) &
        !          + A3(i,2) * B(2,j) &
        !          + A3(i,3) * B(3,j) &
        !          + A3(i,4) * B(4,j) &
        !          + A3(i,5) * B(5,j)
        !C1(i,j) = sum(A1(i,:)*B(:,j))
        !C2(i,j) = sum(A2(i,:)*B(:,j))
        !C3(i,j) = sum(A3(i,:)*B(:,j))
        C1(i,j) = dot_product(A1(i,:),B(:,j))
        C2(i,j) = dot_product(A2(i,:),B(:,j))
        C3(i,j) = dot_product(A3(i,:),B(:,j))
      enddo
    enddo

end subroutine mxm5_3comp

program main
    implicit none
    interface 
    subroutine matmulz(A1,A2,A3,B,C1,C2,C3) bind(c,name="matmulx")
      use iso_c_binding
      implicit none
  
      integer,parameter :: n1 = 25,n3 = 5
      real(kind=c_float),dimension(n1,5),intent(in) :: A1,A2,A3
      real(kind=c_float),dimension(5,n3),intent(in) :: B
      real(kind=c_float),dimension(n1,n3),intent(out) :: C1,C2,C3
      
    end subroutine matmulz 
    end interface 
    integer,parameter :: NGLL =  5
    real :: start, finish
    real(kind=4),dimension(NGLL*NGLL,5) :: A1,A2,A3
    real(kind=4),dimension(5,NGLL) :: B
    real(kind=4),dimension(NGLL*NGLL,NGLL) :: C1,C2,C3
    integer :: i 

    call random_number(B)
    
    call cpu_time(start)
    do i=1,998000
        call random_number(A1)
        call random_number(A2)
        call random_number(A3)
        call mxm5_3comp_singleB(A1,A2,A3,NGLL*NGLL,B,C1,C2,C3,NGLL)
    enddo
    call cpu_time(finish)
    print*,'time elapsed = ',finish - start 

    call cpu_time(start)
    do i=1,998000
        call random_number(A1)
        call random_number(A2)
        call random_number(A3)
        call mxm5_3comp(A1,A2,A3,B,C1,C2,C3)
    enddo
    call cpu_time(finish)
    print*,'time elapsed = ',finish - start 
  
    ! call cpu_time(start)
    ! do i=1,990000
    !     call random_number(A1)
    !     call random_number(A2)
    !     call random_number(A3)
    !     call matmulz(A1,A2,A3,B,C1,C2,C3)
    ! enddo
    ! call cpu_time(finish)
    ! print*,'time elapsed = ',finish - start 

end program 