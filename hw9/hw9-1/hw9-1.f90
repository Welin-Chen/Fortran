Program hw9_1
  implicit none
  integer lon1,lat1,i
  integer pgopen,istat
  real  lon2,lat2
  real  lon(49917),lat(49917),depth(49917),mag(49917)
  real  lon_taiwan(1483),lat_taiwan(1483)
  real  depth_max,mag_max
  REAL  elat,elon,slat,slon,dx,dy,delta
  REAL  b_dx,b_dy,b_delta
  REAL  a_dx,a_dy,a_delta
  integer :: n=1483
  REAL :: x0, y0, x(1483), y(1483),p(49917),d_normal(49917)
  INTEGER :: l, m

  write(*,*)"------------------start------------------"
  Open(1,file='1999.LIS',status='old')
  Open(2,file='output-1999.txt',status='unknown')
if(1==2)then
  open(3,file='Taiwan.txt',status='old')
  Open(4,file='output-Taiwan.txt',status='unknown')
endif

  do i=1,49917
    Read(1,'(18x,i2,f5.2,i3,f5.2,f6.2,f4.2)',end=99)lat1,lat2,lon1,lon2,depth(i),mag(i)
    lon(i)=lon1+lon2/60.0
    lat(i)=lat1+lat2/60.0
    write(2,'(f16.2,2x,f16.2,2x,f16.2,2x,f16.2)')lon(i),lat(i),depth(i),mag(i)
!   write(*,*)lon1,lon2,lat1,lat2,depth,mag
  enddo
  99 close(1)

!以(120,20)為原點,建立B向量
  elat=23.
  elon=120.
  slat=25.
  slon=122.

  call delaz(elat,elon,slat,slon,b_dx,b_dy,b_delta)
  write(*,*)"B vector=",b_dx,b_dy,b_delta

!計算A向量
  do i=1,49917
    call  delaz(elat,elon,lat(i),lon(i),a_dx,a_dy,a_delta)
    write(*,*)"A vector=",i,a_dx,a_dy,a_delta
!p是Ａ在Ｂ的投影量
    p(i)=(a_dx*b_dx+a_dy*b_dy)/sqrt(b_dx**2+b_dy**2)
    d_normal(i)=sqrt(sqrt(a_dx**2+a_dy**2)**2-p(i)**2)
  end do
!繪圖設定
  istat=PGOPEN('hw9-1.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'    call pgslw(2)
  call pgenv(0.0,b_delta,-300.0,0.0,1,1) !just,iaxis
  call pglab('(km)','depth(km)','earthequake distribution')
  call pgsch(2.0)
  do i=1,49917
    if(d_normal(i)<50 .and. p(i)>=0.0 .and. p(i)<=b_delta)then
      call pgsci(2)
      call pgpt1(p(i),-depth(i),1)
    endif
  enddo
  write(*,*)b_delta,maxval(p)

  call pgend()

  write(*,*)"------------------end------------------"


End Program
