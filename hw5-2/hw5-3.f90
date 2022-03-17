Program hw5_3
  implicit none
  integer i,count
  integer pgopen,istat
  real  Distance(5000),Depth(5000),Vp1(5000),Vp_perturbation(5000)
  real  zmin,zmax,xdis,xdep,vp,dx,ivpvs
  real  z,x,y,dfz
  real  xr,xg,xb
  Distance(:)=0
  Depth(:)=0
  Vp1(:)=0
  Vp_perturbation(:)=0


  !馬賽克
  Open(2,file='/Users/chenweilin/desktop/fortran/hw5-2/Vp_prof.dat',status='old')

  count=0
  read(2,*)
  do i=1,5000
    Read(2,'(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)',end=99)Distance(i),Depth(i),Vp1(i),Vp_perturbation(i)
!    write(*,*)Distance(i),Depth(i),Vp(i),Vp_perturbation(i)
    count=count+1
  enddo
99 close(1)
  write(*,*)1,Distance(1),Depth(i),Vp1(1),Vp_perturbation(1)
  write(*,*)count,Distance(count),Depth(count),Vp1(count),Vp_perturbation(count)
  write(*,*) maxval(Depth),maxloc(Depth)

  !繪圖
  istat=PGOPEN('hw5-3.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'

  call pgsubp(1,2)
  !圖1
  call pgslw(4)
  call pgenv(0.0,200.0,-200.0,0.0,1,1) !just,iaxis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Distance(km)','Depth(km)',' vp(km/s)')

  do i=1,4003
    zmin=4.0
    zmax=8.0
    xdis=Distance(i)
    xdep=Depth(i)
    vp=Vp1(i)
    dx=1.0
    call plot_vpvs_value(zmin,zmax,xdis,xdep,vp,dx,ivpvs)
  enddo

  x=100.0
  y=-180.0
  Call pgtext(100.,-190.,'4')
  Call pgtext(130.,-190.,'6')
  Call pgtext(160.,-190.,'8')
  do i=1,60
    zmin=100
    zmax=160
    z=x
    x=x+1
    dfz=(zmax-zmin)/3.0
    if(z.gt.zmax)z=zmax
    if(z.lt.zmin)z=zmin
    if(z.le.(zmin+dfz))then
      xr=1.0
      xg=(z-zmin)/dfz
      xb=0.0
    else if(z.le.(zmin+2.*dfz))then
       xr=1.-(z-zmin-dfz)/dfz
       xg=1.0
       xb=(z-zmin-dfz)/dfz
    else
       xr=0.0
       xg=1.-(z-zmin-2.*dfz)/dfz
       xb=1.0
    endif
    call pgscr(30,xr,xg,xb)
    call pgslw(1)
    call pgsci(30)
    call pgrect(x-dx,x+dx,y-dx,y+dx)
  enddo

  !圖2
  call pgsci(1)
  call pgslw(4)
  call pgenv(0.0,200.0,-200.0,0.0,1,1) !just,iaxis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Distance(km)','Depth(km)',' Vp_perturbation(km/s)')

  ivpvs=1
  do i=1,4003
    zmin=-0.1
    zmax=0.1
    xdis=Distance(i)
    xdep=Depth(i)
    vp=Vp_perturbation(i)*0.01
    dx=1.0
    call plot_vpvs_value(zmin,zmax,xdis,xdep,vp,dx,ivpvs)
  enddo
if(1==2)then
  x=100.0
  y=-180.0
  Call pgtext(100.,-190.,'-10%')
  Call pgtext(160.,-190.,'+10%')
  do i=1,60
    zmin=100
    zmax=160
    z=x
    x=x+1
    dfz=(zmax-zmin)/2.0
    if(z.gt.zmax)z=zmax
    if(z.lt.zmin)z=zmin
    dfz=(zmax-zmin)/2.0
    if(z.le.(zmin+dfz))then
       xr=1.0
       xg=(z-zmin)/dfz
       xb=(z-zmin)/dfz
    else
       xr=1.-(z-zmin-dfz)/dfz
       xg=1.-(z-zmin-dfz)/dfz
       xb=1.0
    endif
    call pgscr(30,xr,xg,xb)
    call pgslw(1)
    call pgsci(30)
    call pgrect(x-dx,x+dx,y-dx,y+dx)
  enddo
endif
  call pgend()
End Program

subroutine plot_vpvs_value(zmin,zmax,xdis,xdep,vp,dx,ivpvs)
  implicit none
  real  zmin,zmax,xdis,xdep,vp,dx,ivpvs
  real  z,x,y,dfz
  real  xr,xg,xb
  z=vp
  x=xdis
  y=-1.*xdep

  dfz=(zmax-zmin)/3.0
  if(z.gt.zmax)z=zmax
  if(z.lt.zmin)z=zmin
  if(ivpvs.eq.1)then
    !-- for vpvs
    dfz=(zmax-zmin)/2.0
    if(z.le.(zmin+dfz))then
       xr=1.0
       xg=(z-zmin)/dfz
       xb=(z-zmin)/dfz
    else
       xr=(zmax-z)/dfz
       xg=(zmax-z)/dfz
       xb=1.0
    endif
  else
!-- for Vp
    if(z.le.(zmin+dfz))then
      xr=1.0
      xg=(z-zmin)/dfz
      xb=0.0
    else if(z.le.(zmin+2.*dfz))then
       xr=1.-(z-zmin-dfz)/dfz
       xg=1.0
       xb=(z-zmin-dfz)/dfz
    else
       xr=0.0
       xg=1.-(z-zmin-2.*dfz)/dfz
       xb=1.0
       if(xg<0.0)then
         write(*,*)xg,(z-zmin-2.*dfz)/dfz
       endif
    endif
  endif
  if(abs(xr)<1e-5)then
    xr=0.0
  endif
  if(abs(xg)<1e-5)then
    xg=0.0
  endif
  if(abs(xb)<1e-5)then
    xb=0.0
  endif

  call pgscr(30,xr,xg,xb)
  call pgslw(1)
  call pgsci(30)
  call pgrect(x-dx,x+dx,y-dx,y+dx)
end subroutine plot_vpvs_value
