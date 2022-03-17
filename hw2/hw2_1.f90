Program hw2_1
  implicit none
  integer a1,a2,a3
  integer b1,b2,b3
  integer Asum,Bsum
  b1=2
  b2=3
  b3=4

  write(*,*)  "猜數字"
  do while(1>0)
    Asum=0
    Bsum=0
10    write(*,*)  "請輸入三個數字a1,a2,a3"
    read(*,*) a1,a2,a3
    write(*,*)  "你輸入的三個數字為",a1,a2,a3
    if(a1==a2 .or. a2==a3 .or. a1==a3)then
      write(*,*)"你輸入的三個數字必須各不相同,請重新輸入"
      goto 10
    endif
    if(a1==b1)then
      Asum=Asum+1
    elseif(a1==b2 .or. a1==b3)then
      Bsum=bsum+1
    endif
    if(a2==b2)then
      Asum=Asum+1
    elseif(a2==b1 .or. a2==b3)then
        Bsum=bsum+1
    endif
    if(a3==b3)then
      Asum=Asum+1
    elseif(a3==b1 .or. a3==b2)then
        Bsum=bsum+1
    endif

    write(*,*) Asum,"A",Bsum,"B"
    if(Asum==3)then
      write(*,*) "恭喜猜對了！"
      exit
    endif
  enddo
End Program
