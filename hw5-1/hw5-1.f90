Program hw5_1
  implicit none
  integer,parameter :: n=30,n1=500
  integer i,j,k
  real,parameter :: V=6.5
  real  x,y,z
  real  xi(n),yi(n),zi(n),di(n),ti(n)
  real  x0,y0,z0
  real  dx,dy,dz
  real  M(4),G(30,4),GT(4,30),GTG(4,4),GTD(4),D(30)
  real  secp(100)
  character(len=10) sta(100),sta_total(n1)
  real  lon(n1),lon1(n1),lon2(n1),lon3(n1),lat(n1),lat1(n1),lat2(n1),lat3(n1),alt(n1),alt3(n1)
  real  elat,elon,slat,slon,delta
  integer len1,len2
  real  t0,dt

  Open(1,file='/Users/chenweilin/Desktop/fortran/hw5-1/21032239.P01',status='old')
  read(1,*)
  len1=0
  do i=1,100
!    Read(1,'(1x,a5,f5.1,13x,f5.1)',end=99)sta(i),secp(i)
    Read(1,'(1x,a4,18x,f6.2)',end=99)sta(i),ti(i)
!    write(*,*)sta(i),ti(i)
    len1=len1+1
  enddo
  99 close(1)
  write(*,*)"len1=",len1

  Open(2,file='/Users/chenweilin/Desktop/fortran/hw5-1/nsta.dat',status='old')
  len2=0
  do i=1,n1
    Read(2,'(a4,f2.0,f5.2,1x,f3.0,f5.2,1x,f6.1)',end=199)sta_total(i),lat1(i),lat2(i),lon1(i),lon2(i),alt(i)
  !  write(*,*)sta_total(i),lat1(i),lat2(i),lon1(i),lon2(i),alt(i)
    lon3(i)=lon1(i)+lon2(i)/60.0
    lat3(i)=lat1(i)+lat2(i)/60.0
    alt3(i)=alt(i)/1000.0
!    write(*,*)sta_total(i),lat(i),lon(i),alt(i)
    len2=len2+1
  enddo
  199 close(2)
  write(*,*)"len2=",len2

  do i=1,len1
    do j=1,len2
      if(sta(i).eq.sta_total(j))then
        lon(i)=lon3(j)
        lat(i)=lat3(j)
        alt(i)=alt3(j)
        exit
      endif
    enddo
!    write(*,*)sta(i),lat(i),lon(i),alt(i)
  enddo

  do i=1,len1
    elat=24.0
    elon=121.0
    slat=lat(i)
    slon=lon(i)
    call delaz(elat,elon,slat,slon,dx,dy,delta)
    xi(i)=dx
    yi(i)=dy
    zi(i)=alt(i)
!    write(*,*) xi(i),yi(i),zi(i),ti(i)
  enddo

  dx=0.
  dy=0.
  dz=0.
  dt=0.
  x0=0.
  y0=0.
  z0=0.
  t0=0.
  !  write(*,*) x0,y0,z0
  write(*,*)V
  write(*,*)"--------start compute--------"
  do k=1,1e4
    x0=x0+dx
    y0=y0+dy
    z0=z0+dz
    t0=t0+dt
    if(z0>0)THEN
      z0=-z0
    endif
    write(*,*)k,"x0=",x0,"y0=",y0,"z0=",z0,"t0=",t0

    do i=1,30
      G(i,1)=xi(i)-x0
      G(i,2)=yi(i)-y0
      G(i,3)=zi(i)-z0
      G(i,4)=-(V**2)*(ti(i)-t0)

      D(i)=0.5*( (xi(i)-x0)**2+(yi(i)-y0)**2+(zi(i)-z0)**2-(V**2)*((ti(i)-t0)**2) )
  !    write(*,*) i,D(i)
    enddo
    GT=transpose(G)

    GTG=matmul(GT,G)
    call MATRIXINV(GTG,4)
    GTD=matmul(GT,D)
    M=matmul(GTG,GTD)

    dx=M(1)
    dy=M(2)
    dz=M(3)
    dt=M(4)
  !    write(*,*)"M", M
    if(abs(dx)<=1e-3 .and. abs(dy)<=1e-3 .and. abs(dz)<=1e-3 .and. abs(dt)<=1e-3)then
      exit
    endif
  enddo
  write(*,*)k,"x0=",x0,"y0=",y0,"z0=",z0,"t0=",t0

!km換成deg
  x0=x0/100.0+121.0
  y0=y0/105.0+24.0
  write(*,*) "lon=",x0,"deg,","lat=",y0,"deg,","alt=",z0,"km"
write(*,*)"--------end compute--------"
End Program
