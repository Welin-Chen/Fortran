real*8 function MakeMStime(iyr,month,iday,ihr,min,second)
! Return input time as MStime--the number of seconds since
! 01/01/1970 00:00:00.0
  integer iyr,month,iday    !Note that year must include century
  integer ihr,min
  real second
  integer basedate    !Julian day for 01/01/1970
  real SperDay          !seconds in a 24 hour day
  integer julianDay					   
  data basedate/2440588/
  data SperDay/86400.0/

  julianDay=JULDAY(month,iday,iyr)	! get astronomical Julian day
  DeltaDays=julianDay-basedate
  MakeMStime=DeltaDays*SperDay+ihr*3600+min*60+second
  return
end function MakeMStime

FUNCTION JULDAY(MM,ID,IYYY)                                  
  PARAMETER (IGREG=15+31*(10+12*1582))                         
  IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'                
  IF (IYYY.LT.0) IYYY=IYYY+1                                   
  IF (MM.GT.2) THEN                                            
    JY=IYYY                                                    
    JM=MM+1                                                    
  ELSE                                                         
    JY=IYYY-1                                                  
    JM=MM+13                                                   
  ENDIF                                                        
  JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995             
  IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN                        
    JA=INT(0.01*JY)                                            
    JULDAY=JULDAY+2-JA+INT(0.25*JA)                            
  ENDIF                                                        
  RETURN                                                       
END function JULDAY