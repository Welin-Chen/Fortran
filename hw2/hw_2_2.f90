Program hw2_2
  implicit none
  integer year,m,hr,min,lon1,lat1,day
  real  sec,lon2,lat2
  Open(1,file='/Users/chenweilin/desktop/fortran/hw2/1999.LIS',status='old')
  Open(2,file='/Users/chenweilin/desktop/fortran/hw2/output.txt',status='unknown')

10   Read(1,'(i4,4i2,f6.2,i2,f5.2,i3,f5.2)',end=99)year,m,day,hr,min,sec,lon1,&
          &lon2,lat1,lat2
     write(2,'(f16.2,2x,f16.2)')lon1+lon2/60.0,lat1+lat2/60.0
     write(*,*)lon1,lon2,lat1,lat2
     goto 10
99 close(1)
End Program
