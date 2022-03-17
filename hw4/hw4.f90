program hw4
  integer pgopen
  real  t(14200),u(14200),n(14200),e(14200)
  real  umax
  integer umaxloc(1),nmaxloc(1),emaxloc(1)
  istat=PGOPEN('test.ps/vcps')  !PostScript
  !-- /c color /v vertical
  if(istat<=0)stop 'ERR opening for PS file!'
!  call pgslw(2)
!  call pgenv(0.0,6.0,0.0,8.0,1,2) !just,iaxis
!  call pgmove(0.0,0.0)
!  call pgdraw(6.0,8.0)
!  call pgend()
  open(1,file='/Users/chenweilin/Desktop/fortran/hw4/seisdata.txt',status='old')

  read(1,*)
  do i=1,14200
    Read(1,'(4F10.3)',end=99)t(i),u(i),n(i),e(i)
  !  write(*,*)t(i),u(i),n(i),e(i)
  enddo
  99 close(1)


  umaxloc=maxloc(abs(u))
  nmaxloc=maxloc(abs(n))
  emaxloc=maxloc(abs(e))
  write(*,*) umaxloc,nmaxloc,emaxloc

  call pgsubp(1,3)
  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,72.0,-100.0,100.0,0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Time (sec)','u','Plot Testing')
  call pgslw(2)
  call pgsci(2)
  call pgline(14200,t,u)
  call pgtext(t(umaxloc),u(umaxloc),'peak')
  call pgsci(6)
  call pgmove(t(umaxloc),-100.)
  call pgdraw(t(umaxloc),100.)

  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,72.0,-100.0,100.0,0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Time (sec)','n','Plot Testing')
  call pgslw(2)
  call pgsci(3)
  call pgline(14200,t,n)
  call pgtext(t(nmaxloc),n(nmaxloc),'peak')
  call pgsci(6)
  call pgmove(t(nmaxloc),-100.)
  call pgdraw(t(nmaxloc),100.)

  call pgsci(1)
  call pgslw(4)
  call pgsch(2.0)
  call pgscf(3)
  call pgenv(0.0,72.0,-200.0,200.0,0,1) !just,axis
  ! Label the axes (note use of \u and \d for raising exponent).
  call pglab('Time (sec)','e','Plot Testing')
  call pgslw(2)
  call pgsci(4)
  call pgline(14200,t,e)
  call pgtext(t(emaxloc),e(emaxloc),'peak')
  call pgsci(6)
  call pgmove(t(emaxloc),-200.)
  call pgdraw(t(emaxloc),200.)
  call pgend()

end program hw4
