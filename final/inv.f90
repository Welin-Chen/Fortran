program main
  implicit none
  real*8  A(3,3),inv_A(3,3)
  data A / 1,2,3,4,5,6,7,8,8/
  call inverse_matrix(A,inv_A)
end

subroutine inverse_matrix(A,inv_A)
  implicit none
  integer i,j
  real*8  A(3,3),inv_A(3,3),I_matrix(3,3)
  real*8  f
  data I_matrix / 1,0,0,0,1,0,0,0,1 /

!  write(*,*)  "A="
!  write(*,200)((A(i,j),j=1,3),i=1,3)
  200 format(3f12.2)

  !計算反矩陣
  !1
  if(A(1,1) /= 1)then
    f=1.0/A(1,1)
    do j=1,3
      A(1,j)=f*A(1,j)
      I_matrix(1,j)=f*I_matrix(1,j)
    enddo
  endif
  !2
  if(A(2,1) /= 0)then
    f=-A(2,1)/A(1,1)
    do j=1,3
      A(2,j)=A(2,j)+f*A(1,j)
      I_matrix(2,j)=I_matrix(2,j)+f*I_matrix(1,j)
    enddo
  endif
  !3
  if(A(3,1) /= 0)then
    f=-A(3,1)/A(1,1)
    do j=1,3
      A(3,j)=A(3,j)+f*A(1,j)
      I_matrix(3,j)=I_matrix(3,j)+f*I_matrix(1,j)
    enddo
  endif
  !4
  if(A(2,2) /= 1)then
    f=1.0/A(2,2)
    do j=1,3
      A(2,j)=f*A(2,j)
      I_matrix(2,j)=f*I_matrix(2,j)
    enddo
  endif
  !5
  if(A(3,2) /= 0)then
    f=-A(3,2)/A(2,2)
    do j=1,3
      A(3,j)=A(3,j)+f*A(2,j)
      I_matrix(3,j)=I_matrix(3,j)+f*I_matrix(2,j)
    enddo
  endif
  !6
  if(A(3,3) /= 1)then
    f=1.0/A(3,3)
    do j=1,3
      A(3,j)=f*A(3,j)
      I_matrix(3,j)=f*I_matrix(3,j)
    enddo
  endif
  !7
  if(A(2,3) /= 0)then
    f=-A(2,3)/A(3,3)
    do j=1,3
      A(2,j)=A(2,j)+f*A(3,j)
      I_matrix(2,j)=I_matrix(2,j)+f*I_matrix(3,j)
    enddo
  endif
  !8
  if(A(1,3) /= 0)then
    f=-A(1,3)/A(3,3)
    do j=1,3
      A(1,j)=A(1,j)+f*A(3,j)
      I_matrix(1,j)=I_matrix(1,j)+f*I_matrix(3,j)
    enddo
  endif
  !9
  if(A(1,2) /= 0)then
    f=-A(1,2)/A(2,2)
    do j=1,3
      A(1,j)=A(1,j)+f*A(2,j)
      I_matrix(1,j)=I_matrix(1,j)+f*I_matrix(2,j)
    enddo
  endif
  inv_A=I_matrix
  do i=1,3
    do j=1,3
      if(abs(inv_A(i,j))>1e10)then
        write(*,*)"inverse matrix dosen't exit!"
        stop
      endif
    enddo
  enddo
!  write(*,*)  "inv_A="
!  write(*,200)((inv_A(i,j),j=1,3),i=1,3)
  return
end
