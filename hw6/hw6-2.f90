Program hw6_2
  implicit none
  real  dist(100),az(100),toa(100),min_1(100),secp(100),resp(100),&
        &wtp(100),secs(100),ress(100),wts(100),smag(100)
  real  secp_new(100)
  real  s,sx,sy,sxx,sxy,delta,a,b,r
  real  sum1,sum21,sum22,sum3,var,sigma
  real  xavg,yavg,yi
  character sta(100)
  integer i
  integer pgopen,istat
  real  x1,y1,x0,y0
  Character card*20
  real  Distance(5000),Depth(5000),Vp(5000),Vp_perturbation(5000)
  integer ndata
  real rms,r_cor

  Open(1,file='/Users/chenweilin/desktop/fortran/hw6/ppfile.txt',status='old')

  read(1,*)
  do i=1,100
  Read(1,'(6x,f5.1,13x,f5.1)',end=99)dist(i),secp(i)
  enddo

99 close(1)
  sx=0
  sy=0
  sxx=0
  sxy=0
  delta=0
  a=0
  b=0
  do i=1,30
    secp_new(i)=secp(i)-69.23
    write(*,*)dist(i),secp(i),secp_new(i)
    s=s+1
    sx=sx+dist(i)
    sy=sy+secp_new(i)
    sxx=sxx+dist(i)**2
    sxy=sxy+dist(i)*secp_new(i)
  enddo

  delta=s*sxx-sx**2
  a=(sxx*sy-sx*sxy)/delta
  b=(s*sxy-sx*sy)/delta
  xavg=sx/30.0
  yavg=sy/30.0
  do i=1,30
    sum1=sum1+(dist(i)-xavg)*(secp_new(i)-yavg)
  !  sum21=sum21+sqrt((dist(i)-xavg)**2)*sqrt((secp_new(i)-yavg)**2)
    sum21=sum21+(dist(i)-xavg)**2
    sum22=sum22+(secp_new(i)-yavg)**2
  end do
  r=sum1/(sqrt(sum21)*sqrt(sum22))

  sum3=0
  do i=1,30
    yi=a+b*dist(i)
    sum3=sum3+(secp_new(i)-yi)**2
  end do
  var=sum3/29
  sigma=sqrt(var)

  write(*,*)  "s=",s
  write(*,*)  "sx=",sx
  write(*,*)  "sy=",sy
  write(*,*)  "sxx=",sxx
  write(*,*)  "sxy=",sxy
  write(*,*)  "delta=",delta
  write(*,*)  "a=",a
  write(*,*)  "b=",b
  write(*,*)  "r=",r
  write(*,*)  "var=",var
  write(*,*)  "sigma=",sigma


!    call xyfit(x,y,ndata,cept,slop,rms,r_cor,sdv)

  ndata=30
  call xyfit(dist,secp_new,ndata,a,b,rms,r_cor,sigma)
  !繪圖
  istat=PGOPEN('hw6-2.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'
  call pgslw(2)
  call pgenv(0.0,120.0,0.0,20.0,1,1) !just,iaxis
!    call pgmove(0.0,0.0)
!    call pgdraw(6.0,8.0)
  call pglab('dist(km)','secp_new(sec)','earthequake')
  call pgsci(2)



  x0=0
  y0=a
  x1=120
  y1=a+b*x1
  call pgsci(1)
  call pgmove(x0,y0)
  call pgdraw(x1,y1)

  call pgsls(2)
  call pgsci(3)
  call pgmove(x0,y0+sigma)
  call pgdraw(x1,y1+sigma)

  call pgsls(2)
  call pgsci(4)
  call pgmove(x0,y0-sigma)
  call pgdraw(x1,y1-sigma)

  do i=1,30
    call pgsci(2)
    call pgpt1(dist(i),secp_new(i),20)
  enddo

!  call pgsci(2)
!  call pgtext(119.,25.2,'<depth_max*(1/3)')

  Write(card,10)a,b
10 format('Y = ',f6.2,' + ',f6.2,'X')
  Call pgtext(70.,5.,trim(card))

  Write(card,20)r_cor
20 format('R = ',f6.2)
  Call pgtext(10.,10.,trim(card))

  Write(card,30)sigma
30 format('sdv = ',f6.2)
  Call pgtext(10.,15.,trim(card))

End Program
