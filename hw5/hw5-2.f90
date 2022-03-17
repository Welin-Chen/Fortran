Program hw5
  implicit none
  integer lon1,lat1,i
  integer pgopen,istat
  real  lon2,lat2
  real  lon(49917),lat(49917),depth(49917),mag(49917)
  real  lon_taiwan(1483),lat_taiwan(1483)
  real  depth_max,mag_max

  Open(1,file='/Users/chenweilin/Desktop/fortran/hw5/1999.LIS',status='old')
  Open(2,file='/Users/chenweilin/desktop/fortran/hw5/output-1999.txt',status='unknown')
  open(3,file='/Users/chenweilin/Desktop/fortran/hw5/Taiwan.txt',status='old')
  Open(4,file='/Users/chenweilin/desktop/fortran/hw5/output-Taiwan.txt',status='unknown')

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

  istat=PGOPEN('hw5.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'
    call pgslw(2)
    call pgenv(119.0,123.0,21.0,26.0,1,1) !just,iaxis
!    call pgmove(0.0,0.0)
!    call pgdraw(6.0,8.0)
    call pglab('longitude(E)','latitude(N)','earthequake distribution')
    call pgsci(2)
    call pgtext(119.,25.2,'<depth_max*(1/3)')
    call pgsci(3)
    call pgtext(119.,25.4,'<depth_max*(2/3)')
    call pgsci(4)
    call pgtext(119.,25.6,'>depth_max*(2/3)')

    do i=1,49917
      mag_max=maxval(mag)
      depth_max=maxval(depth)
    enddo
    write(*,*)"mag_max=",mag_max
    write(*,*)"depth_max=",depth_max
    do i=1,49917
call pgsci(3)
      if(depth(i)<depth_max/3.0)then
        call pgsci(2)
      elseif(depth(i)<depth_max/3.0*2.0)then
        call pgsci(3)
      else
        call pgsci(4)
      endif
      if(mag(i)<mag_max/3.0)then
        call pgpt1(lon(i),lat(i),20)
      elseif(mag(i)<mag_max/3.0*2.0)then
        call pgpt1(lon(i),lat(i),22)
      else
        call pgpt1(lon(i),lat(i),24)
      endif
    enddo
    call pgsci(1)
    do i=1,1483
      call pgpt1(lon_taiwan(i),lat_taiwan(i),1)
    enddo
    call pgend()
End Program
