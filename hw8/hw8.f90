Program hw8
  implicit none
  integer i,j
  integer num
  integer pgopen,istat
  character A(832),ax(3),ay(3),az(3)
  real  axd(4600),ayd(4600),azd(4600)
  real  dt,t(4600)
  real  sum,avg
  character*23 cdate

  type :: header
    character*4 code
    character*4 dummy_char
    real*8 orgintime
    integer*2 ncom
    integer*2 dummy_int
    integer*4 ndata
    real dt
    real dummy_real
  end type header
  type(header)::wh
  real,allocatable :: wd(:,:)

  write(*,*)"----------start----------"

  Open(1,file='/Users/chenweilin/Desktop/fortran/hw8/seismicdata.bin',form='unformatted',status='old',access="stream")

  read(1)wh
  allocate(wd(wh%ncom,wh%ndata))
  write(*,*)wh%ncom,wh%ndata
  read(1)wd

!  write(*,*)wh%orgintime,wh%ncom,wh%ndata
  !write(*,*)wd
  do i=1,wh%ndata
    axd(i)=wd(1,i)
    ayd(i)=wd(2,i)
    azd(i)=wd(3,i)
    t(i)=i*wh%dt
    write(*,*)t(i),axd(i),ayd(i),azd(i)
  enddo
  call decode_mstime(cdate,wh%orgintime)
  write(*,*)cdate
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
  istat=PGOPEN('hw8.ps/vcps')  !PostScript
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
