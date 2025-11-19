MODULE base_module
  IMPLICIT NONE
  ! common parameter
  REAL(8)::lx,ly                 ! domain size
  INTEGER::nxc,nyc               ! number of cells
  REAL(8)::dx,dy                 ! grid spacing

  contains 
    FUNCTION xc(i) result (x)
      INTEGER::i
      REAL(8)::x
      x = (0.5d0+dble(i-1))*dx
    END FUNCTION xc
    FUNCTION yc(j) result (y)
      INTEGER::j
      REAL(8)::y
      y = (0.5d0+dble(j-1))*dy - 0.5*ly
    END FUNCTION yc
    FUNCTION tact() result(tpr)
      integer :: seed, itmp
      real :: tpr,rnd,area
      real :: c1, mm  ! NSD parameters
      ! Set the seed value
      seed = 1234
      c1 = 2.549d7
      mm = 1.1451d0
      area = lx * ly
      call random_seed(seed)
      call random_number(rnd)
      itmp = 1 + int(rnd * ((nxc * nyc)-1))
      tpr = (itmp/(area*c1))**(1.0/mm)
    END FUNCTION tact
END MODULE

PROGRAM SITE
  USE base_module
  IMPLICIT NONE
  INTEGER::i,j,gLevel
  REAL(8)::tact_sum,tmp,tact_max,tact_min

  gLevel = 32  !32
  lx = 0.01d0
  ly = 0.005d0
  nxc = 16*gLevel
  nyc = 8*gLevel
  dx = lx/dble(nxc)
  dy = ly/dble(nyc)
  tact_sum=0.0d0
  tact_min=1.0d10
  tact_max=0.0d10

  !WRITE(*,*)xc(nxc),lx-xc(nxc),xc(1)
  !WRITE(*,*)yc(nyc),ly-yc(nyc),yc(1)
  OPEN(10,file='site.txt')
    WRITE(10,*)"X      Y      Tact"
    DO j = 1,nyc
    DO i = 1,nxc
      tmp=tact()
      WRITE(10,*)xc(i),yc(j),tmp
      tact_sum = tact_sum + tmp
      IF (tact_min > tmp) tact_min = tmp
      IF (tact_max < tmp) tact_max= tmp
    ENDDO
    ENDDO
  CLOSE (10)
  WRITE(*,*)"nxc, nyc=",nxc, nyc
  WRITE(*,*)"dx, dy=",dx, dy
  WRITE(*,*)"Tact min max=",tact_min, tact_max
  WRITE(*,*)"Total number of nucleation site =",nxc*nyc
  WRITE(*,'(a,E12.3)')" nucleation site density in defined =",dble(nxc*nyc)/(lx*ly)
  WRITE(*,*)"Average of Tact=",tact_sum/dble(nxc*nyc)

  STOP
END PROGRAM SITE

