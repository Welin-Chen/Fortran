Program hw5_3
  implicit none
  integer i,j
  integer num
  integer pgopen,istat
  character A(832),ax(3),ay(3),az(3)
  real  axd(4600),ayd(4600),azd(4600)
  real  dt,t(4600)
  real  sum,avg

  write(*,*)"----------start----------"

  Open(1,file='/Users/chenweilin/Desktop/fortran/hw5-3/03134001.CVA',form='unformatted',status='old',access="stream")
  do i=1,832
    read(1)A(i)
  enddo
  write(*,*)Ichar(A(33)),Ichar(A(34)),Ichar(A(35)),Ichar(A(36))
  num=Ichar(A(33))*256**3+Ichar(A(34))*256**2+Ichar(A(35))*256**1+Ichar(A(36))
  write(*,*)"num=",num
  write(*,*)Ichar(A(41)),Ichar(A(42))
  dt=Ichar(A(41))*256+Ichar(A(42))
  dt=dt*0.001
  write(*,*)"dt=",dt,"sec"

  do i=1,4600
    do j=1,3
      read(1)ax(j)
    enddo
    do j=1,3
      read(1)ay(j)
    enddo
    do j=1,3
      read(1)az(j)
    enddo
    call bitodec(ax,ay,az,axd,ayd,azd,i)
    t(i)=i*dt
    write(*,*)t(i),axd(i),ayd(i),azd(i)
  enddo
  sum=0.0
  do i=1,100
    sum=sum+axd(i)
  enddo
  avg=sum/100.0
  do i=1,4600
    axd(i)=axd(i)-avg
  enddo

  sum=0.0
  do i=1,100
    sum=sum+ayd(i)
  enddo
  avg=sum/100.0
  do i=1,4600
    ayd(i)=ayd(i)-avg
  enddo

  sum=0.0
  do i=1,100
    sum=sum+azd(i)
  enddo
  avg=sum/100.0
  do i=1,4600
    azd(i)=azd(i)-avg
  enddo
!繪圖
  istat=PGOPEN('hw5-3.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'
  call pgsubp(1,3)
  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,50.0,-10.0,10.0,0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Time (sec)','n(cm/s^2)','Plot Testing')
  call pgslw(2)
  call pgsci(2)
  call pgline(4600,t,axd)

  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,50.0,-10.0,10.0,0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Time (sec)','e(cm/s^2)','Plot Testing')
  call pgslw(2)
  call pgsci(2)
  call pgline(4600,t,ayd)

  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,50.0,-10.0,10.0,0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Time (sec)','u(cm/s^2)','Plot Testing')
  call pgslw(2)
  call pgsci(2)
  call pgline(4600,t,azd)

  call pgend()

  write(*,*)"----------end----------"


End Program
subroutine bitodec(ax,ay,az,axd,ayd,azd,i)
  integer i
  character ax(3),ay(3),az(3)
  real  axd(4600),ayd(4600),azd(4600)
  axd(i)=Ichar(ax(1))*256**2+Ichar(ax(2))*256+Ichar(ax(1))
  ayd(i)=Ichar(ay(1))*256**2+Ichar(ay(2))*256+Ichar(ay(1))
  azd(i)=Ichar(az(1))*256**2+Ichar(az(2))*256+Ichar(az(1))
  if(axd(i)>2**23)then
    axd(i)=axd(i)-2**24
  endif
  if(ayd(i)>2**23)then
    ayd(i)=ayd(i)-2**24
  endif
  if(azd(i)>2**23)then
    azd(i)=azd(i)-2**24
  endif
  axd(i)=axd(i)*2*980./(2**23)
  ayd(i)=ayd(i)*2*980./(2**23)
  azd(i)=azd(i)*2*980./(2**23)
end
