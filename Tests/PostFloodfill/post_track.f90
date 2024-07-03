  PROGRAM POST_TRACK
  IMPLICIT NONE

!------------------------------------------------------------------------------!
  TYPE :: Bubble
    INTEGER :: bid, bcvol  ! bubble ID and cvol
    REAL(8) :: bx, by, bz  ! bubble x, y, z
    INTEGER :: iactive     ! 1 for active, 0 for deactive
    INTEGER :: ndata       ! number of time step data
    REAL(8) :: cvol_ave    ! cvol averaged in time
  END TYPE Bubble
!------------------------------------------------------------------------------!
  TYPE(Bubble), ALLOCATABLE :: bubbles(:)
  INTEGER :: nb_total      ! total number of bubbles
  INTEGER :: nb_current    ! number of bubble at current time step
  INTEGER,ALLOCATABLE :: cvol(:),ID(:),i2ib(:)  ! array for current time step
  REAL(8),ALLOCATABLE :: x(:),y(:),z(:)         ! array for current time step
  REAL(8) :: tol_dist      ! Tolerance for distance. If difference of bubble
                           !  center positions between previous and current time
                           !  step is smaller than tol_dist, then the bubble is
                           !  assumed to be same.
  INTEGER :: tol_cvol      ! Tolerance for cvol. Bubbles below tol_cvol will be
                           !  neglected.
  REAL(8) :: tol_zmax      ! Tolerance for max z. If bubble z-position is below
                           !  tol_zmax, the bubble will be neglected in 
                           !  statistics and gnuplot
  INTEGER :: tol_step      ! Tolerance for time steps. If a bubble exists
                           !  shorter than tol_step (number of time step),
                           !  then the bubble will be neglected in statistics
                           !  and gnuplot
!------------------------------------------------------------------------------!
  ! temporary
  INTEGER :: i, j, ib        ! counter
  INTEGER :: idtmp,icvol
  REAL(8) :: t_in,t_old    
  REAL(8) :: ttmp,xyz(3),dist
  CHARACTER(len=64) :: ctmp
  CHARACTER(len=512) :: fmt
  ! for gnuplot
  INTEGER :: nb_plot
  INTEGER,ALLOCATABLE :: ib_plot(:)
!------------------------------------------------------------------------------!
  tol_dist = 1e-3
  tol_cvol = 300
  tol_zmax = 0.02
  tol_step = 100

  ! initialize ----------------------------------------------------------------!
  nb_total = 0
  ! crude code: max number of bubbles are limited here
  ALLOCATE(bubbles(10000))
  bubbles%ndata = 0;
  bubbles%cvol_ave = 0.0;

  ! Preprocess ----------------------------------------------------------------!
  CALL system("grep Info1 tracked_regions.txt > info1.out")
  WRITE(*,*)"Output data to info1.out"
  WRITE(*,*)"Output data to info1_negID.out"
  OPEN(10,file="info1.out")
  OPEN( 9,file="info1_negID.out")
    DO i = 1, 999999
      READ(10,*,err=990,end=991)ctmp,ttmp,ctmp,idtmp,ctmp,icvol,ctmp,(xyz(j),j=1,3)
      IF (idtmp<0.and.icvol>tol_cvol) THEN
        WRITE(9,'(i5,E16.8,i7,3E16.8)')idtmp,ttmp,icvol,(xyz(j),j=1,3)
      ENDIF
    ENDDO
    990 WRITE(*,*)"File read error:info1.out line=",i
    991 CONTINUE
  CLOSE(10)
  CLOSE(9)

  ! Main process --------------------------------------------------------------!
  OPEN(10,file="info1_negID.out")

    500 CONTINUE
    !WRITE(*,*)"Pass 500"

    READ(10,*,end=999,err=999)idtmp,t_old
    nb_current = 1
    DO i = 1, 10000
      READ(10,*,end=998,err=998)idtmp,t_in
      !IF (dabs(t_old-0.351666)<1e-7) THEN
      ! WRITE(*,*)t_in,idtmp
      !ENDIF
      IF (t_old.eq.t_in) THEN
      !IF (dabs(t_old-t_in)<1e-7) THEN
        nb_current = nb_current +1
      ELSE
        !WRITE(*,*)"HELLO",nb_current
        EXIT
      ENDIF
    ENDDO

    998 CONTINUE

    DO i = 1, nb_current+1
      BACKSPACE (10)
    ENDDO

    ALLOCATE(x(nb_current),y(nb_current),z(nb_current))
    ALLOCATE(ID(nb_current),cvol(nb_current),i2ib(nb_current))

    DO i = 1, nb_current
      READ(10,*)ID(i),t_in,cvol(i),x(i),y(i),z(i)
    ENDDO

    IF (nb_total.eq.0) THEN

      ! first step
      bubbles(1)%bid = 1
      bubbles(1)%bcvol = cvol(1)
      bubbles(1)%bx = x(1)
      bubbles(1)%by = y(1) 
      bubbles(1)%bz = z(1) 
      bubbles(1)%iactive = 1
      nb_total = 1

    ELSE

      WRITE(*,*)nb_current,t_in
      ! check if same bubble?
      i2ib=0
      DO i = 1, nb_current
        DO ib = 1, nb_total
          IF (bubbles(ib)%iactive==1) THEN
          dist = dsqrt((bubbles(ib)%bx - x(i))**2  &
                      +(bubbles(ib)%by - y(i))**2  &
                      +(bubbles(ib)%bz - z(i))**2)
          IF(dist < tol_dist)THEN
            i2ib(i) = ib
            !WRITE(*,*)i2ib(i)
          ELSE
            !WRITE(*,*)"dist= ",dist,nb_current
          ENDIF
        ENDIF
        ENDDO
      ENDDO

      ! reset iactive
      DO ib = 1, nb_total
        bubbles(ib)%iactive=0
      ENDDO

      ! copy
      DO i = 1, nb_current
        IF (i2ib(i).eq.0) THEN
          nb_total = nb_total+1
          bubbles(nb_total)%bid = nb_total
          bubbles(nb_total)%bcvol = cvol(i)
          bubbles(nb_total)%bx = x(i)
          bubbles(nb_total)%by = y(i)
          bubbles(nb_total)%bz = z(i)
          bubbles(nb_total)%iactive = 1
        ELSE
          bubbles(i2ib(i))%bcvol = cvol(i)
          bubbles(i2ib(i))%bx = x(i)
          bubbles(i2ib(i))%by = y(i)
          bubbles(i2ib(i))%bz = z(i)
          bubbles(i2ib(i))%iactive = 1
        ENDIF
      ENDDO
    ENDIF

    DEALLOCATE(x,y,z,cvol,ID,i2ib)

    !check
    DO ib = 1, nb_total
      IF (bubbles(ib)%iactive.eq.1) THEN
        WRITE(100+ib,'(E16.8,i4,i6,3E16.8)')t_in,bubbles(ib)%bid,  &
          bubbles(ib)%bcvol,bubbles(ib)%bx,bubbles(ib)%by,bubbles(ib)%bz
        bubbles(ib)%ndata = bubbles(ib)%ndata+1
      ENDIF
      IF (bubbles(ib)%bid==0) STOP
    ENDDO

    GOTO 500

    999 CONTINUE

  CLOSE(10)

  ! statistics ----------------------------------------------------------------!
  OPEN(10,file="cvol_ave.out")
  OPEN(11,file="cvol_ave_all.out")
  WRITE(*,*)"Output time-averaged cvol to 'cvol_ave.out' and 'cvol_ave_all.out'"
  DO ib = 1, nb_total
    OPEN(100+ib)
      REWIND(100+ib)
      DO i = 1, bubbles(ib)%ndata
        READ(100+ib, *)ttmp,idtmp,icvol
        bubbles(ib)%cvol_ave = bubbles(ib)%cvol_ave + icvol
      ENDDO
    CLOSE(100+ib)
    bubbles(ib)%cvol_ave = bubbles(ib)%cvol_ave / bubbles(ib)%ndata
    IF (bubbles(ib)%bz>tol_zmax .and. bubbles(ib)%ndata>tol_step) THEN
      WRITE(*,*)"ID=",ib,"ndata= ",bubbles(ib)%ndata,"cvol_ave= ",bubbles(ib)%cvol_ave
      WRITE(10,*)"ID=",ib,"ndata= ",bubbles(ib)%ndata,"cvol_ave= ",bubbles(ib)%cvol_ave
    ENDIF
    WRITE(11,*)"ID=",ib,"ndata= ",bubbles(ib)%ndata,"cvol_ave= ",bubbles(ib)%cvol_ave
  ENDDO
  CLOSE(10)

  ! gnuplot -------------------------------------------------------------------!
  ! set nb_plot
  nb_plot = 0
  DO ib = 1, nb_total
    IF (bubbles(ib)%bz>tol_zmax .and. bubbles(ib)%ndata>tol_step) THEN
      nb_plot = nb_plot+1
    ENDIF
  ENDDO
  WRITE(*,*)'nb_plot=',nb_plot

  ALLOCATE(ib_plot(nb_plot))
  i = 1
  DO ib = 1, nb_total
    IF (bubbles(ib)%bz>tol_zmax .and. bubbles(ib)%ndata>tol_step) THEN
      ib_plot(i)=ib
      i = i+1
      WRITE(*,*)'ib= ',ib,'i=',i
    ENDIF
  ENDDO

  WRITE(*,*)"Output bubble Z position to 'bubbles-z.png'"
  OPEN(10,file='bubbles-z.gnu')
  WRITE(10,*)"set term png size 1600,1200"
  WRITE(10,*)'set output "bubbles-z.png"'
  WRITE(10,*)"set gri"
  WRITE(10,*)'set xlabel "Time (s)"'
  WRITE(10,*)'set ylabel "Z (m)"'
  WRITE(10,*)'l=6'
  WRITE(fmt, '(A,I3,A)') '(a,', nb_plot, '(a,i3,a))'
  !WRITE(*,*)'fmt= ',fmt
  WRITE(10,fmt)'plot ',('"fort.',100+ib_plot(ib),'" u 1:l w l,',ib=1,nb_plot)
  WRITE(10,*)
  Close(10)

  CALL system("gnuplot bubbles-z.gnu")

  STOP
  END
