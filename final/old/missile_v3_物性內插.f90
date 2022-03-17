!------------------------------------------------------------!
!  6Dof彈道模擬(6Degreed of freedom trajectory simulation)    !
!  受控體(plant):無控火箭                                      !
!  作者:陳維霖 台大機械博士班 學號:D09522006                     !
!  date:2021/1/16                                            !
!  version:1                                                 !
!  輸入檔:(1)input_1_thrust                                   !
!        (2)input_2_physical                                   !
!        (3)input_3_aero                                     !
!        (3-1-1)input_3_aero_1_1_CA_ON                       !
!        (3-1-2)input_3_aero_1_2_CA_OFF                      !
!        (3-1-3)input_3_aero_1_3_CN                          !
!        (3-1-4)input_3_aero_1_4_XCP                         !
!        (3-1-5)input_3_aero_1_5_CLFI                        !
!        (3-1-6)input_3_aero_1_6_CHMAL                       !
!        (3-1-7)input_3_aero_1_7_CHMAW                       !
!        (3-1-8)input_3_aero_1_8_CNFAL                       !
!        (3-1-9)input_3_aero_1_9_CNFAW                       !
!        (3-2)input_3_aero_2_reference                       !
!        (3-3)input_3_aero_3_other                           !
!        (4)input_4_atmosphere                               !
!  輸出檔:(1)output_1_thrust_check                            !
!        (2)output_2_physical_check                          !
!        (3)output_3_aero_check                              !
!        (4)output_4_atmosphere_check                        !
!------------------------------------------------------------!
!------------------------------------------------------------!
!                       建模 for 常數、變數                     !
!------------------------------------------------------------!
module constant
  implicit none
  real*8, parameter :: PI=3.1415927410125732421875
  real*8, parameter :: g=9.80665
endmodule

module global
  implicit none
  real*8, save :: t,dt
  integer, save :: counter,i_main,call_times
endmodule

module module_atmosphere
  implicit none
  real*8, save :: P_air,rho,Cs
endmodule

module solution
  implicit none
  real*8, save :: Y(13),dY_dt(13)
  real*8, save :: position(3),velocity(3),angular_velocity(3),Euler_angle_GD(3),mass
  real*8, save :: angular_velocity_rate(3),Euler_angle_rate_GD(3)
  real*8, save :: velocity_D(3)
  real*8, save :: AF,phi_c,apha,V,Ma
endmodule

module input_data
  implicit none
  real*8  input_1_thrust(100,100)
  real*8  input_2_physical(10,10)
  real*8  input_3_aero_1(11,10,9),input_3_aero_2(4),input_3_aero_3(13,9)
endmodule
!------------------------------------------------------------!
!                            主程式                           !
!------------------------------------------------------------!
Program main
  use input_data
  implicit none
  write(*,*)"--------------------------6Dof程式start--------------------------"


  !-------------------讀取輸入資料(input)-------------------!
  !呼叫副程式:讀取輸入資料(input)
  call sub_input()
  !-------------------初始資料計算(initial)-------------------!
  !呼叫副程式:初始資料計算(initial)
  call sub_initial()
  write(*,*)"--------------------------6Dof程式end--------------------------"

End Program
!------------------------------------------------------------!
!  副程式:讀取輸入資料(input)                                   !
!  輸入:                                                      !
!  輸出:input_1_thrust                                        !
!      input_2_physical                                      !
!      input_3_aero_1                                        !
!      input_3_aero_2                                        !
!      input_3_aero_3                                        !
!  定義:input_1_thrust:推力                                   !
!      input_2_physical:物性                                 !
!      input_3_aero_1:氣動力係數(二維)                         !
!      input_3_aero_2:氣動力參考長度,面積                       !
!      input_3_aero_3:氣動力係數(一維)                         !
!      D3:矩陣個數                                            !
!      row:矩陣列數                                           !
!      column:矩陣行數                                        !
!------------------------------------------------------------!
subroutine sub_input()
  use input_data
  integer, parameter :: unit_thrust=10,unit_thrust_check=15
  integer, parameter :: unit_physical=20,unit_physical_check=25
  integer, parameter :: unit_aero_3_1=30,unit_aero_3_2=40,unit_aero_3_3=50
  integer, parameter :: unit_aero_check=140
  integer, parameter :: unit_aero=111
  integer i,j,k
  integer row,column,d3
  integer status

  write(*,*)"------------------副程式:讀取輸入資料(input)------------------"

  !(1)讀取資料-推力
  row=26
  column=5
  Open(unit=unit_thrust,file='input_1_thrust.txt',status='old')
    read(unit_thrust,*)
    read(unit_thrust,*,iostat=status) ((input_1_thrust(i,j),j=1,column),i=1,row)
    write(*,*)"推力資料,讀取結果(0:成功,other:失敗):",status
  close(unit_thrust)

  Open(unit=unit_thrust_check,file='output_1_thrust_check.txt',status='unknown')
    write(unit_thrust_check,"(5A16)") "t(sec)","p(kgf/cm^2)","TF_sea(kgf)","TF_vacuum(kgf)","dm/dt(kg/sec)"
    write(unit_thrust_check,"(5f16.2)") ((input_1_thrust(i,j),j=1,column),i=1,row)
  close(unit_thrust_check)

  !(2)讀取資料-物性
  row=2
  column=4
  Open(unit=unit_physical,file='input_2_physical.txt',status='old')
    read(unit_physical,*)
    read(unit_physical,*,iostat=status) ((input_2_physical(i,j),j=1,column),i=1,row)
    write(*,*)"物性資料,讀取結果(0:成功,other:失敗):",status
  close(unit_physical)

  Open(unit=unit_physical_check,file='output_2_physical_check.txt',status='unknown')
    write(unit_physical_check,"(5A16)") "m(kg)","x_cg(m)","Ixx(kg-m^2)","Iyy(kg-m^2)"
    write(unit_physical_check,"(4f16.2)") ((input_2_physical(i,j),j=1,column),i=1,row)
  close(unit_physical_check)

  !(3)讀取資料-氣動力係數
  row=11
  column=10
  d3=9
  !3-1
  !1-1:CA(ON)
  k=1
  Open(unit=unit_aero_3_1,file='input_3_aero_1_1.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-2:CA(OFF)
  k=2
  Open(unit=unit_aero_3_1,file='input_3_aero_1_2.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-3:CN
  k=3
  Open(unit=unit_aero_3_1,file='input_3_aero_1_3.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-4:XCP
  k=4
  Open(unit=unit_aero_3_1,file='input_3_aero_1_4.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-5:CLFI(22.5)
  k=5
  Open(unit=unit_aero_3_1,file='input_3_aero_1_5.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-6:CHMAL
  k=6
  Open(unit=unit_aero_3_1,file='input_3_aero_1_6.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-7:CHMAW
  k=7
  Open(unit=unit_aero_3_1,file='input_3_aero_1_7.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-8:CNFAL
  k=8
  Open(unit=unit_aero_3_1,file='input_3_aero_1_8.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !1-9:CNFAW
  k=9
  Open(unit=unit_aero_3_1,file='input_3_aero_1_9.txt',status='old')
    read(unit_aero_3_1,*)
    read(unit_aero_3_1,*,iostat=status) ((input_3_aero_1(i,j,k),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_1)

  !3-2
  k=10
  Open(unit=unit_aero_3_2,file='input_3_aero_2.txt',status='old')
    read(unit_aero_3_2,"(9x,f16.8)",iostat=status) (input_3_aero_2(i),i=1,4)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_2)

  !3-3
  row=13
  column=9
  k=11
  Open(unit=unit_aero_3_3,file='input_3_aero_3.txt',status='old')
    read(unit_aero_3_3,*,iostat=status) ((input_3_aero_3(i,j),j=1,column),i=1,row)
    write(*,*)"氣動力係數資料k=",k,",讀取結果(0:成功,other:失敗):",status
  close(unit_aero_3_3)

  !檢查氣動力輸入資料
  Open(unit=unit_aero_check,file='output_3_aero_check.txt',status='unknown')
  !3-1
  row=11
  column=10
  d3=9
  do k=1,9
      if(k==1)then
        write(unit_aero_check,*)k,":CA(ON)"
      elseif(k==2)then
        write(unit_aero_check,*)k,":CA(OFF)"
      elseif(k==3)then
        write(unit_aero_check,*)k,":CN"
      elseif(k==4)then
        write(unit_aero_check,*)k,":XCP"
      elseif(k==5)then
        write(unit_aero_check,*)k,":CLFI"
      elseif(k==6)then
        write(unit_aero_check,*)k,":CHMAL"
      elseif(k==7)then
        write(unit_aero_check,*)k,":CHMAW"
      elseif(k==8)then
        write(unit_aero_check,*)k,":CNFAL"
      elseif(k==9)then
        write(unit_aero_check,*)k,":CNFAW"
      endif
      write(unit_aero_check,100) ("MACH=",(input_3_aero_1(i,j,k),j=1,column),i=1,1)
      write(unit_aero_check,100) ("AF=",(input_3_aero_1(i,j,k),j=1,column),i=2,row)
      write(unit_aero_check,"(/)")
  enddo
  !3-2
  write(unit_aero_check,101) "RL(M)=",input_3_aero_2(1)
  write(unit_aero_check,101) "RA(M^2)=",input_3_aero_2(2)
  write(unit_aero_check,101) "ALT(KM)=",input_3_aero_2(3)
  write(unit_aero_check,"(a9,e16.4)") "RE=",input_3_aero_2(4)
  write(unit_aero_check,"(/)")
  !3-3
  row=13
  column=9
!  write(unit_aero_check,"(9f10.2)") ((input_3_aero_3(i,j),j=1,column),i=1,row)
  i=1
  write(unit_aero_check,102) i,"MACH=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CNQ=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CMQ=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CLP=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CND=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"XCPD=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CLDP=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CHMD=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CNFD=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"XCPFD=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CND2=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"XCPD2=",(input_3_aero_3(i,j),j=1,column)
  i=i+1
  write(unit_aero_check,102) i,"CLDP2=",(input_3_aero_3(i,j),j=1,column)

  100 format(a8,10f8.2)
  101 format(a9,f16.8)
  102 format(i3,a8,9f8.2)

  close(unit_aero_check)
  write(*,*)"----------------------------------------------------------"
  return
end
!------------------------------------------------------------!
!  副程式:初始資料計算(initial)                                 !
!  輸入:input_1_thrust                                       !
!      input_2_physical                                      !
!  輸出:position                                             !
!      thrust                                               !
!      physical                                             !
!      dt                                                   !
!      count                                                !
!  定義:position(1)=x=飛彈在地面座標系g之x軸位置(km)             !
!      position(2)=y=飛彈在地面座標系g之y軸位置(km)             !
!      position(3)=z=飛彈在地面座標系g之z軸位置(km)             !
!      velocity(1)=u=飛彈在地面座標系g之x軸速度(m/s)            !
!      velocity(2)=v=飛彈在地面座標系g之y軸速度(m/s)            !
!      velocity(3)=w=飛彈在地面座標系g之z軸速度(m/s)            !
!      angular_velocity(1)=p=飛彈在地面座標系g之x角速度(rad/s)  !
!      angular_velocity(2)=q=飛彈在地面座標系g之y角速度(rad/s)  !
!      angular_velocity(3)=r=飛彈在地面座標系g之z角速度(rad/s)  !
!      Euler_angle_GD(1)=theta=俯仰角(rad)                   !
!      Euler_angle_GD(2)=phi=偏航角(rad)                     !
!      Euler_angle_GD(3)=gamma=傾斜角(rad)                   !
!      thrust=內插完的推力                                    !
!      physical=內插完的物性                                  !
!      dt=時間步長(s)                                        !
!      counter=到燃畢時的計算次數                              !
!------------------------------------------------------------!
subroutine sub_initial()
  use global
  use constant
  use input_data
  use solution
  implicit none
  integer, parameter :: unit_thrust_check=15
  integer, parameter :: unit_physical_check=25
  integer i,j,k
  integer row,column
  real*8 thrust_test(5000,5)
  real*8 physical_test(5000,4)
  real*8 t_start,t_end
  real*8 A_up,A_mid,A_down
  real*8 B_up,B_mid,B_down
  real*8 mass_up,mass_down,mass_loss,mass_loss_rate
  real*8 mass_start,mass_end
  real*8 factor
  real*8 L_body,X_CG
  real*8 t_up,t_down

  write(*,*)"----------------副程式:初始資料計算(initial)----------------"
  dt=0.01
  !L_body=彈體長度(m)
  L_body=3.873
  !X_CG=彈體重心位置(m)
  X_CG=2.69

  !初始俯仰角=25度
  Euler_angle_GD(1)=25.0*PI/180.
  Euler_angle_GD(2)=0.0
  Euler_angle_GD(3)=0.0

  position(1)=0.0
  position(2)=0.0
  position(3)=-0.8

  velocity(1)=0.0
  velocity(2)=0.0
  velocity(3)=0.0

  angular_velocity(1)=0.0
  angular_velocity(2)=0.0
  angular_velocity(3)=0.0

  write(*,"(a4,f10.4,a4)")"dt=",dt,"(s)"
  write(*,300)"L_body(m)","X_CG(m)","L_body-X_CG(m)"
  write(*,310)L_body,X_CG,L_body-X_CG
  write(*,300)"theta(deg)","phi(deg)","gamma(deg)"
  write(*,310)(Euler_angle_GD(i),i=1,3)
  write(*,300)"x(m)","y(m)","z(m)"
  write(*,310)(position(i),i=1,3)
  write(*,300)"Vx(m/s)","Vy(m/s)","Vz(m/s)"
  write(*,310)(velocity(i),i=1,3)
  write(*,300)"p(deg/s)","q(deg/s)","r(deg/s)"
  write(*,310)(angular_velocity(i)*PI/180.,i=1,3)

  300 format(3A20)
  310 format(3F20.4)
  !--------------------內插--------------------!
  !1.推力內插
  counter=(input_1_thrust(26,1)-input_1_thrust(1,1))/dt+1
  write(*,*)"內插筆數:",counter
  t_start=input_1_thrust(1,1)
  t_end=input_1_thrust(26,1)
  row=counter
  column=5
  do j=1,5
    do i=1,row
      t=t_start+(i-1)*dt
      if((t-t_start)<1e-5)then
        thrust_test(i,j)=input_1_thrust(1,j)
      elseif(abs(t-t_end)<=1e-5)then
        thrust_test(i,j)=input_1_thrust(26,j)
      else
        do k=2,26
          t_up=input_1_thrust(k-1,1)
          t_down=input_1_thrust(k,1)
          if(t>=t_up .and. t<=t_down)then
            A_up   = t_up
            A_mid  = t
            A_down = t_down
            B_up   = input_1_thrust(k-1,j)
            B_down = input_1_thrust(k,j)
            B_mid  = B_up + ( (A_mid-A_up)/(A_down-A_up) )*(B_down-B_up)
            thrust_test(i,j) = B_mid
          endif
        enddo
      endif
    end do
  end do
  !輸出確認
  open(unit=unit_thrust_check,file="output_thrust_內插.txt")
    write(unit_thrust_check,100) "t(sec)","p(kgf/cm^2)","TF_sea(kgf)","TF_vacuum(kgf)","dm/dt(kg/sec)"
    write(unit_thrust_check,110)((thrust_test(i,j),j=1,column),i=1,row)
    100 format(5a16)
    110 format(5f16.2)
  close(unit_thrust_check)

  !2.物性內插
  row = counter
  column = 4

  factor = 1.0
  do while(.true.)
    mass_start = input_2_physical(1,1)
    mass_end = input_2_physical(2,1)
    mass_loss = 0.0
    do i=1,row
      mass_loss_rate     = factor*thrust_test(i,5)
      mass_loss          = mass_loss + mass_loss_rate*dt
      mass               = mass_start - mass_loss
      physical_test(i,1) = mass
    end do
    if( abs(physical_test(row,1)-mass_end) <= 1e-2)then
      exit
    elseif((physical_test(row,1)-mass_end)>0)then
      factor = factor + (physical_test(row,1)-mass_end)/mass_end
    else
      factor = factor + (physical_test(row,1)-mass_end)/mass_end
    endif
  end do
  write(*,300)"mass_start(kg)","mass_end_compute","mass_end"
  write(*,310)mass_start,mass,mass_end
  write(*,"(a40,f20.4)")"mass_end_compute-mass_end(kg)=",physical_test(row,1)-mass_end

  do j=2,4
    physical_test(1,j)   = input_2_physical(1,j)
    physical_test(row,j) = input_2_physical(2,j)
    mass_up   = physical_test(1,1)
    mass_down = physical_test(row,1)
    A_up      = mass_up
    A_down    = mass_down
    B_up      = input_2_physical(1,j)
    B_down    = input_2_physical(2,j)
    do i=1,row-1
      A_mid = physical_test(i,1)
      B_mid = B_up + ( (A_mid-A_up)/(A_down-A_up) )*(B_down-B_up)
      physical_test(i,j) = B_mid
    end do
  end do
  !輸出確認
  open(unit=unit_physical_check,file="output_physical_內插.txt")
    write(unit_physical_check,200) "m(kg)","x_cg(m)","Ixx(kg-m^2)","Iyy(kg-m^2)"
    write(unit_physical_check,210)((physical_test(i,j),j=1,column),i=1,row)
    200 format(4a16)
    210 format(4f16.2)
  close(unit_physical_check)

  !設定初始值
  !積分的變數(dY)共13個:1~3:速度、4~6:位置、7~9:角速度、10~12:尤拉角、13:質量
  do i=1,3
    Y(i)   = velocity(i)
    Y(i+3) = position(i)
    Y(i+6) = angular_velocity(i)
    Y(i+9) = Euler_angle_GD(i)
  enddo
  Y(13) = input_2_physical(1,1)
  mass  = Y(13)
  write(*,"(i16,f16.2)")(i,Y(i),i=1,13)
  write(*,*)"----------------------------------------------------------"
  return
end
