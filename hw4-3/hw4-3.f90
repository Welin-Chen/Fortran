Program hw4_3
  implicit none
  integer,parameter :: n=5
  integer i,j,k
  real  x,y,z
  real  xi(n),yi(n),zi(n),di(n)
  real  x0,y0,z0
  real  dx,dy,dz
  real  M(3),G(5,3),GT(3,5),GTG(3,3),GTD(3),D(5)

  i=1
  xi(i)=1.
  yi(i)=1.
  zi(i)=0.
  di(i)=1.42

  i=2
  xi(i)=-1.
  yi(i)=-1.
  zi(i)=0.
  di(i)=1.42

  i=3
  xi(i)=1.
  yi(i)=-1.
  zi(i)=0.
  di(i)=1.41

  i=4
  xi(i)=-1.
  yi(i)=1.
  zi(i)=0.
  di(i)=1.43

  i=5
  xi(i)=2.
  yi(i)=0.
  zi(i)=0.
  di(i)=2.1

  write(*,*) "i,xi(i),yi(i),zi(i),di(i)"
  do i=1,5
    write(*,*) i,xi(i),yi(i),zi(i),di(i)
  enddo
!  (1,1,0,1.42),(-1,-1,0,1.42),(1,-1,0,1.41),
!  (-1,1,0,1.43),(2,0,0,2.1)

  x0=20.
  y0=20.
  z0=20.
!  write(*,*) x0,y0,z0

  write(*,*)"--------start compute--------"
  do k=1,1e6
    x0=x0+dx
    y0=y0+dy
    z0=z0+dz
      write(*,"(i6,3(a6,e16.4))")k,"x0=",x0,"y0=",y0,"z0=",z0

    do i=1,5
      G(i,1)=x0-xi(i)
      G(i,2)=y0-yi(i)
      G(i,3)=z0-zi(i)

      D(i)=0.5*(di(i)**2-(xi(i)-x0)**2-(yi(i)-y0)**2-(zi(i)-z0)**2)
  !    write(*,*) i,D(i)
    enddo
    GT=transpose(G)

    GTG=matmul(GT,G)
    call MATRIXINV(GTG,3)
    GTD=matmul(GT,D)
    M=matmul(GTG,GTD)

    dx=M(1)
    dy=M(2)
    dz=M(3)
!    write(*,*)"M", M
    if(abs(dx)<=1e-7 .and. abs(dy)<=1e-7 .and. abs(dz)<=1e-7)then
      exit
    endif

  enddo
  write(*,"(i6,3(a6,e16.4))")k,"x0=",x0,"y0=",y0,"z0=",z0
End Program
