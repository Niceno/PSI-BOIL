  Program pathline
  implicit none
  integer::nstart,nend,interval,nstep,istep,itime,ndigit
  integer::np,i
  real(8),allocatable::x(:,:),y(:,:),z(:,:),u(:,:),v(:,:),w(:,:),time(:)
  character(len=512),allocatable::fname(:)
  character(len=512)::fncommon
  character(len=512)::ctmp,ctmp1


  ndigit=6

  101 CONTINUE
  WRITE(*,*)"Input start of time step"
  READ(*,*,err=101,end=101) nstart

  102 CONTINUE
  WRITE(*,*)"Input end of time step"
  READ(*,*,err=102,end=102) nend

  103 continue
  write(*,*)"Input time step increment"
  read(*,*,err=103,end=103)interval

  104 continue
  write(*,*)"Input common file name."
  write(*,*)"Example: particles_"
  read(*,*,err=104,end=104)fncommon

  nstep=1+(nend-nstart)/interval

  ALLOCATE(fname(nstep))

  ! set file names
  DO istep = 1,nstep
    itime = nstart+(istep-1)*interval
    call int2char(itime,ctmp,ndigit)
    fname(istep) = trim(fncommon)//ctmp(1:ndigit)//".dat"
    !WRITE(*,*)trim(fname(istep))
  ENDDO

  ! set np
  OPEN(10,file=trim(fname(1)),status='old',err=999)
    READ(10,*)
    READ(10,*)ctmp,ctmp,np
  CLOSE(10)

  ! allocate variables
  allocate(x(np,nstep),y(np,nstep),z(np,nstep), &
           u(np,nstep),v(np,nstep),w(np,nstep),time(nstep))

  ! read data
  DO istep = 1, nstep
    OPEN(10,file=trim(fname(istep)),status='old',err=998)
    READ(10,*)
    READ(10,*)
    !READ(10,*)ctmp,time(istep)
    DO i = 1, np
      READ(10,*)x(i,istep),y(i,istep),z(i,istep), &
                u(i,istep),v(i,istep),w(i,istep)
    ENDDO
    CLOSE(10)
    WRITE(*,*)trim(fname(istep))
  ENDDO

  OPEN(10,file='pathline.dat',FORM='formatted')
    DO i = 1, np
      call int2char(i,ctmp1,ndigit)
      ctmp = 'TITLE= "Particle'//ctmp1(1:ndigit)//'"'
      !WRITE(10,*)'TITLE= "Particle',i,'"'
      WRITE(10,*)trim(ctmp)
      WRITE(10,*)"VARIABLES= X Y Z U V W"
      WRITE(10,*)"ZONE I=",nstep, "DATAPACKING=POINT"
      DO istep = 1, nstep
        WRITE(10,'(6E16.8)')x(i,istep),y(i,istep),z(i,istep), &
                u(i,istep),v(i,istep),w(i,istep)
      ENDDO
    ENDDO
  CLOSE(10)

  STOP

 998 continue
   WRITE(*,*)"Error!!! File not found. ",trim(fname(istep))
   STOP

 999 continue
   WRITE(*,*)"Error!!! File not found. ",trim(fname(1))
   STOP

  END


!-----------------------------------------------------------------------
   SUBROUTINE int2char( num, cwrk, ndigit )
   INTEGER num, j
   CHARACTER(LEN=1024) cwrk
   CHARACTER(LEN=1) CH1
!
   DO j=1,ndigit
     ch1=char(48+mod(num/10**(ndigit-j),10))
     cwrk(j:j)=ch1
   ENDDO

   RETURN
   END SUBROUTINE int2char

