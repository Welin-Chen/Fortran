Program hw3_1
  implicit none
  real  dist(100),az(100),toa(100),min_1(100),secp(100),resp(100),&
        &wtp(100),secs(100),ress(100),wts(100),smag(100)
  real  secp_new(100)
  real  s,sx,sy,sxx,sxy,delta,a,b,r
  real  sum1,sum21,sum22,sum3,var,sigma
  real  xavg,yavg,yi
  character sta(100)
  integer i

  Open(1,file='/Users/chenweilin/desktop/fortran/hw3/ppfile.txt',status='old')

  read(1,*)
  do i=1,100
    Read(1,'(6x,f5.1,13x,f5.1)',end=99)dist(i),secp(i)
  enddo

  pause

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


End Program
