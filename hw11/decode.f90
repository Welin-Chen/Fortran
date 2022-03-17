!-----------------------------------------------------------
!       SUDS TIME CODE COUNTS SINCE 1970 IN SECOND
!       IN REAL*8 DATA TYPE SAVE. I TRANS IT TO CHAR*23
!-----------------------------------------------------------
subroutine decode_mstime(cdate,long_real)
    character*23 cdate
    real*8 long_real
    idate=idint(long_real/86400.)
    xtime=long_real-idate*86400.
	idate=idate+1
    ih=int(xtime/3600.)
    im=int(xtime/60.)-ih*60
    ss=xtime-im*60.-ih*3600.
    ic=0
    do i=1970,9999
      id=i-int(i/100)*100
      if(id.eq.0)then
        id=int(i/100)
        id=id-int(id/4)*4
        if(id.eq.0)id=100
      else
        id=id-int(id/4)*4
        if(id.eq.0)id=100
      endif
      idd=365
      if(id.eq.100)idd=366
      ic=ic+idd
      if(idate.le.ic)then
        idate=idate-ic+idd
        goto 10
      endif
    enddo
10  iyear=i
    if(idate.lt.31)then
      imo=1
      idat=idate
      goto 20
    endif
    jdate=31
    ifre=28
    if(id.eq.100)ifre=29
    jdate1=jdate+ifre
    if(idate.gt.jdate.and.idate.le.jdate1)then
      imo=2
      idat=idate-jdate
      goto 20
    endif
    do k=3,12
      jdate=jdate1
      if(k.eq.3.or.k.eq.5.or.k.eq.7.or.k.eq.8.or.k.eq.10.or.k.eq.12)ifre=31
      if(k.eq.4.or.k.eq.6.or.k.eq.9.or.k.eq.11)ifre=30
      jdate1=jdate+ifre
      if(idate.gt.jdate.and.idate.le.jdate1)then
        imo=k
        idat=idate-jdate
        goto 20
      endif
    enddo
20  cdate='mm/dd/yyyy hh:mm:ss.sss'

	write(cdate,30)imo,idat,iyear,ih,im,ss
30  format(i2,'/',i2,'/',i4,1x,i2,':',i2,':',f6.3)
    return
end subroutine decode_mstime	