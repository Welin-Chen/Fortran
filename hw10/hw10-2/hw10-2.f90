Program hw9_2
  implicit none
  integer i,j
  integer num
  integer pgopen,istat
  character A(832),ax(3),ay(3),az(3)
  real  axd(4600),ayd(4600),azd(4600)
  real  dt,t(4600)
  real  sum,avg
  REAL  XREAL(4096),XIMAG(4096),axd_mag(4096),ayd_mag(4096),azd_mag(4096),f,fi(4096)
  REAL  axd_mag_2(4096),ayd_mag_2(4096),azd_mag_2(4096)
  real  axd_2(4600),ayd_2(4600),azd_2(4600)
  INTEGER N,NU

  write(*,*)"----------start----------"

  Open(1,file='03134001.CVA',form='unformatted',status='old',access="stream")
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

  dt=0.01
  N=4096
  NU=12
  f=1./dt
  write(*,*)N*dt,f
  do i=1,N
    XREAL(i)=axd(i)
    XIMAG(i)=0.0
  enddo
  !axd_mag
  call FFT(XREAL,XIMAG,N,NU)
  do i=1,N
    if(i==1)then
      fi(i)=1.0/(N*dt)
    else
      fi(i)=1.0/(N*dt)+fi(i-1)
    endif
    axd_mag(i)=sqrt(XREAL(i)**2+XIMAG(i)**2)/4096.0
    write(*,*)i,axd_mag(i)
!0.5~3hz保留
    if( (fi(i)>=0.5 .and. fi(i)<=3.) .or. (fi(i)<=(100.-0.5) .and. fi(i)>=(100.-3.)) )then
      axd_mag_2(i)=XREAL(i)
      XIMAG(i)=-XIMAG(i)
    else
      axd_mag_2(i)=0.0
      XIMAG(i)=0.0
    endif
  enddo
  !轉回時域
  call FFT(axd_mag_2,XIMAG,N,NU)
  do i=1,N
    axd_2(i)=axd_mag_2(i)/4096.0
  end do

  !ayd_mag
  do i=1,N
    XREAL(i)=ayd(i)
    XIMAG(i)=0.0
  enddo
  call FFT(XREAL,XIMAG,N,NU)
  do i=1,N
    ayd_mag(i)=sqrt(XREAL(i)**2+XIMAG(i)**2)/4096.0
  enddo
  !azd_mag
  do i=1,N
    XREAL(i)=azd(i)
    XIMAG(i)=0.0
  enddo
  call FFT(XREAL,XIMAG,N,NU)
  do i=1,N
    azd_mag(i)=sqrt(XREAL(i)**2+XIMAG(i)**2)/4096.0
  enddo

!繪圖
  istat=PGOPEN('hw9-2.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'
  call pgsubp(1,3)
  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
!  call pgenv(0.0,50.0,-10.0,10.0,0,1) !just,axis
  call pgenv(0.0,50.0,0.0,MAXVAL(axd_mag),0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('f (hz)','axd_mag','Plot Testing')
  call pgslw(2)
  call pgsci(2)
!  call pgline(4600,t,axd)
  call pgline(4096/2,fi,axd_mag)

  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,50.0,0.0,MAXVAL(ayd_mag),0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('f (hz)','ayd_mag','Plot Testing')
  call pgslw(2)
  call pgsci(2)
  call pgline(4096/2,fi,ayd_mag)


  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,50.0,0.0,MAXVAL(azd_mag),0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('f (hz)','azd_mag','Plot Testing')
  call pgslw(2)
  call pgsci(2)
  call pgline(4096/2,fi,azd_mag)


  call pgend()
  !繪圖
    istat=PGOPEN('hw10-2.ps/vcps')  !PostScript
    !-- /c color /v vertical
    if(istat<=0)stop 'ERR opening for PS file!'
    call pgsubp(1,4)
    call pgsci(1)
    call pgslw(4)
    call pgsch(2.0)
    call pgscf(3)
  !  call pgenv(0.0,50.0,-10.0,10.0,0,1) !just,axis
    call pgenv(0.0,50.0,0.0,MAXVAL(axd_mag),0,1) !just,axis
    ! Label the axes (note use of \u and \d for raising exponent).
    call pglab('f (hz)','axd_mag','Plot Testing')
    call pgslw(2)
    call pgsci(2)
  !  call pgline(4600,t,axd)
    call pgline(4096/2,fi,axd_mag)

    call pgsci(1)
    call pgslw(4)
    call pgsch(2.0)
    call pgscf(3)
  !  call pgenv(0.0,50.0,-10.0,10.0,0,1) !just,axis
    call pgenv(0.0,50.0,0.0,MAXVAL(axd_mag_2),0,1) !just,axis
    ! Label the axes (note use of \u and \d for raising exponent).
    call pglab('f (hz)','axd_mag_2','Plot Testing')
    call pgslw(2)
    call pgsci(2)
  !  call pgline(4600,t,axd)
    call pgline(4096/2,fi,axd_mag_2)

    call pgsci(1)
    call pgslw(4)
    call pgsch(2.0)
    call pgscf(3)
    call pgenv(0.0,50.0,-MAXVAL(axd),MAXVAL(axd),0,1) !just,axis
    ! Label the axes (note use of \u and \d for raising exponent).
    call pglab('t (sec)','axd_mag(filter)','Plot Testing')
    call pgslw(2)
    call pgsci(2)
    call pgline(4096,t,axd)


    call pgsci(1)
    call pgslw(4)
    call pgsch(2.0)
    call pgscf(3)
    call pgenv(0.0,50.0,-MAXVAL(axd_2),MAXVAL(axd_2),0,1) !just,axis
    ! Label the axes (note use of \u and \d for raising exponent).
    call pglab('t (sec)','axd_mag(filter)','Plot Testing')
    call pgslw(2)
    call pgsci(2)
    call pgline(4096,t,axd_2)

  !  call FFT(XREAL,XIMAG,N,NU)
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
