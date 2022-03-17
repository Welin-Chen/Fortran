program main
  implicit none
  integer i,j
  real  ma(3,3),mi(3,3),inv_ma(3,3)
  real  f
!  data ma / 2,1,2,3,4,6,6,9,1/
  data ma / 2,1,4,-1,0,0,3,-2,2/
  data mi / 1,0,0,0,1,0,0,0,1 /
  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do

write(*,*)  "--------------1--------------"
if(ma(1,1) /= 1)then
  f=1.0/ma(1,1)
  write(*,*)  "f=",f

  do j=1,3
    ma(1,j)=f*ma(1,j)
    mi(1,j)=f*mi(1,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------2--------------"
if(ma(2,1) /= 0)then
  f=-ma(2,1)/ma(1,1)
  write(*,*)  "f=",f
  do j=1,3
    ma(2,j)=ma(2,j)+f*ma(1,j)
    mi(2,j)=mi(2,j)+f*mi(1,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------3--------------"
if(ma(3,1) /= 0)then
  f=-ma(3,1)/ma(1,1)
  write(*,*)  "f=",f
  do j=1,3
    ma(3,j)=ma(3,j)+f*ma(1,j)
    mi(3,j)=mi(3,j)+f*mi(1,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------4--------------"
if(ma(2,2) /= 1)then
  f=1.0/ma(2,2)
  write(*,*)  "f=",f

  do j=1,3
    ma(2,j)=f*ma(2,j)
    mi(2,j)=f*mi(2,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------5--------------"
if(ma(3,2) /= 0)then
  f=-ma(3,2)/ma(2,2)
  write(*,*)  "f=",f
  do j=1,3
    ma(3,j)=ma(3,j)+f*ma(2,j)
    mi(3,j)=mi(3,j)+f*mi(2,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------6--------------"
if(ma(3,3) /= 1)then
  f=1.0/ma(3,3)
  write(*,*)  "f=",f

  do j=1,3
    ma(3,j)=f*ma(3,j)
    mi(3,j)=f*mi(3,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------7--------------"
if(ma(2,3) /= 0)then
  f=-ma(2,3)/ma(3,3)
  write(*,*)  "f=",f
  do j=1,3
    ma(2,j)=ma(2,j)+f*ma(3,j)
    mi(2,j)=mi(2,j)+f*mi(3,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------8--------------"
if(ma(1,3) /= 0)then
  f=-ma(1,3)/ma(3,3)
  write(*,*)  "f=",f
  do j=1,3
    ma(1,j)=ma(1,j)+f*ma(3,j)
    mi(1,j)=mi(1,j)+f*mi(3,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
      write(*,*)  mi(i,1),mi(i,2),mi(i,3)
  end do
endif
write(*,*)  "--------------9--------------"
if(ma(1,2) /= 0)then
  f=-ma(1,2)/ma(2,2)
  write(*,*)  "f=",f
  do j=1,3
    ma(1,j)=ma(1,j)+f*ma(2,j)
    mi(1,j)=mi(1,j)+f*mi(2,j)
  enddo

  write(*,*)  "ma="
  do i=1,3
      write(*,*)  ma(i,1),ma(i,2),ma(i,3)
  end do
  write(*,*)  "inv_ma="
  do i=1,3
    !  write(*,*)  mi(i,1),mi(i,2),mi(i,3)
    write(*,*)  mi(i,j=(1,3))
  end do
endif
end program
