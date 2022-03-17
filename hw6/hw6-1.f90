Program hw6_1
  implicit none
  integer lon1,lat1,i
  integer pgopen,istat
  real  lon2,lat2
  real  lon(49917),lat(49917),depth(49917),mag(49917)
  real  lon_taiwan(1483),lat_taiwan(1483)
  real  depth_max,mag_max
  integer :: n=1483
  REAL :: x0, y0, x(1483), y(1483)
  INTEGER :: l, m

  Open(1,file='1999.LIS',status='old')
  Open(2,file='output-1999.txt',status='unknown')
  open(3,file='Taiwan.txt',status='old')
  Open(4,file='output-Taiwan.txt',status='unknown')

  do i=1,49917
    Read(1,'(18x,i2,f5.2,i3,f5.2,f6.2,f4.2)',end=99)lat1,lat2,lon1,lon2,depth(i),mag(i)
    lon(i)=lon1+lon2/60.0
    lat(i)=lat1+lat2/60.0
    write(2,'(f16.2,2x,f16.2,2x,f16.2,2x,f16.2)')lon(i),lat(i),depth(i),mag(i)
!   write(*,*)lon1,lon2,lat1,lat2,depth,mag
  enddo
  99 close(1)

  do i=1,1483
    Read(3,'(f7.3,f7.3)',end=199)lon_taiwan(i),lat_taiwan(i)
    write(4,'(f16.3,2x,f16.3)')lon_taiwan(i),lat_taiwan(i)
  enddo
  199 close(1)

  do i=1,1483
    x(i)=lon_taiwan(i)
    y(i)=lat_taiwan(i)
  enddo
  n=1483

  istat=PGOPEN('hw6.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'
    call pgslw(2)
    call pgenv(119.0,123.0,21.0,26.0,1,1) !just,iaxis
!    call pgmove(0.0,0.0)
!    call pgdraw(6.0,8.0)
    call pglab('longitude(E)','latitude(N)','earthequake distribution')
    call pgsci(2)
    call pgtext(119.,25.2,'')
    call pgsci(1)
    do i=1,1483
      call pgpt1(lon_taiwan(i),lat_taiwan(i),1)
    enddo
    call pgsci(2)
    do i=1,49917
      x0=lon(i)
      y0=lat(i)
      call locpt (x0,y0,x,y,n,l,m)
      if(l==1)then
    !    write(*,*)i,x0,y0,l
        call pgpt1(x0,y0,1)
      endif
    enddo
    call pgend()




End Program
