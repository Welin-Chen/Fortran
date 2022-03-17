Program hw10_1
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
  REAL  p1(3),p2(3)
  REAL  p12_dx,p12_dy,p12_dz,p12_delta
  REAL  n12_dx,n12_dy
  REAL  dip_dx,dip_dy,dip_dz
  REAL  p12_dx_n,p12_dy_n,p12_dz_n
  REAL  dip_dx_n,dip_dy_n,dip_dz_n
  REAL  p_dx,p_dy,p_delta
  REAL  p_s(49917),p_dip(49917),d(49917)
  INTEGER sum
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
    lon(i)=real(lon1)+lon2/60.0
    lat(i)=real(lat1)+lat2/60.0
    write(2,'(f16.2,2x,f16.2,2x,f16.2,2x,f16.2)')lon(i),lat(i),depth(i),mag(i)
!   write(*,*)lon1,lon2,lat1,lat2,depth,mag
  enddo
  99 close(1)

!以(120,20)為原點,建立B向量
  elat=23.6
  elon=120.6
  slat=24.3
  slon=120.8
!  elat=23.
!  elon=120.
!  slat=25.
!  slon=122.

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
  if(istat<=0)stop 'ERR opening for PS file!'
  call pgslw(2)
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

!hw10-2算在斷層面的投影點及距離小於10km
  p1(1)=120.6
  p1(2)=23.6
  p1(3)=0.
  p2(1)=120.8
  p2(2)=24.3
  p2(3)=0.
!以p1(120.6,23.6)為原點,p2(120.8,24.3)建立p12向量,即走向
  elat=23.6
  elon=120.6
  slat=24.3
  slon=120.8

  call delaz(elat,elon,slat,slon,p12_dx,p12_dy,p12_delta)
  p12_dz=0.0
  write(*,*)"--------------------"
  write(*,*)"p12 vector=",p12_dx,p12_dy,p12_dz
  write(*,*)"p12 vector length=",sqrt(p12_dx**2+p12_dy**2+p12_dz**2)
  write(*,*)"--------------------"
  p12_dx_n=p12_dx/sqrt(p12_dx**2+p12_dy**2)
  p12_dy_n=p12_dy/sqrt(p12_dx**2+p12_dy**2)
  p12_dz_n=0.0
  write(*,*)"p12 vector n=",p12_dx_n,p12_dy_n,p12_dz_n
  write(*,*)"p12 vector n length=",sqrt(p12_dx_n**2+p12_dy_n**2+p12_dz_n**2)
  write(*,*)"--------------------"

!與走向p12垂直(內積為零)且平行水平面的向量
  n12_dx=p12_dy_n
  n12_dy=-p12_dx_n
!  n12_dx=p12_dy
!  n12_dy=-p12_dx
  write(*,*)"n12 vector n=",n12_dx,n12_dy
  write(*,*)"--------------------"

!傾向dip:垂直走向且落於斷層面的向量,斷層面dip長40km,該向量與水平面傾角40度
!  dip_dx=n12_dx*(40.*cos(40.*3.14159/180.0)/p12_delta)
!  dip_dy=n12_dy*(40.*cos(40.*3.14159/180.0)/p12_delta)
  dip_dx=n12_dx*40.*cos(40.*3.14159/180.0)
  dip_dy=n12_dy*40.*cos(40.*3.14159/180.0)
  dip_dz=40.*sin(40.*3.14159/180.0)
  write(*,*)"dip vector=",dip_dx,dip_dy,dip_dz
  write(*,*)"dip vector length=",sqrt(dip_dx**2+dip_dy**2+dip_dz**2)
  write(*,*)"--------------------"

  dip_dx_n=dip_dx/sqrt(dip_dx**2+dip_dy**2+dip_dz**2)
  dip_dy_n=dip_dy/sqrt(dip_dx**2+dip_dy**2+dip_dz**2)
  dip_dz_n=dip_dz/sqrt(dip_dx**2+dip_dy**2+dip_dz**2)
  write(*,*)"dip vector n=",dip_dx_n,dip_dy_n,dip_dz_n
  write(*,*)"dip vector length=",sqrt(dip_dx_n**2+dip_dy_n**2+dip_dz_n**2)
  write(*,*)"--------------------"

  write(*,*)"dip_n與strike_n是否垂直?兩者內積=",dip_dx_n*p12_dx_n+dip_dy_n*p12_dy_n
  write(*,*)"--------------------"

    do i=1,49917
      !計算p向量

      call  delaz(elat,elon,lat(i),lon(i),p_dx,p_dy,p_delta)
!      write(*,*)elat,elon
    !  write(*,*)"p vector=",i,p_dx,p_dy,p_delta
      !p_s:p在strike走向的投影長度
      p_s(i)=(p_dx*p12_dx_n+p_dy*p12_dy_n)
      !p_dip:p在dip傾向的投影長度
      p_dip(i)=(p_dx*dip_dx_n+p_dy*dip_dy_n+dip_dz_n*depth(i))
      !d:p點與斷層面之垂直距離
      d(i)=sqrt((p_dx**2+p_dy**2+depth(i)**2)-(p_s(i)**2+p_dip(i)**2))
    end do

    !繪圖設定
      istat=PGOPEN('hw10-1.ps/vcps')  !PostScript
      !-- /c color /v vertical
      if(istat<=0)stop 'ERR opening for PS file!'
      call pgslw(2)
      call pgenv(0.0,p12_delta,40.0,0.0,1,1) !just,iaxis
      call pglab('(km)','depth(km)','earthequake distribution')
      call pgsch(2.0)
      sum=0.0
      do i=1,49917
        if(d(i)<=10.0)then
          call pgsci(2)
          call pgpt1(p_s(i),p_dip(i),1)
          sum=1+sum
        endif
      enddo
      write(*,*)sum

      call pgend()

  write(*,*)"------------------end------------------"


End Program
