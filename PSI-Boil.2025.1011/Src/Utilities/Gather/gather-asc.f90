!  gather tecplot (PSI-Boil dat format) file
   module base_var
     IMPLICIT NONE
     REAL(4),ALLOCATABLE::x(:,:,:,:),v(:,:,:,:)
     INTEGER::icmax,jcmax,kcmax   !!cell
     INTEGER::inmax,jnmax,knmax   !!node
     INTEGER::np,nvariable,ndigit
     CHARACTER(len=2048),allocatable::fname(:)
     CHARACTER(len=2048)fncommon
   end module base_var

!-----------------------------------------------------------------------
 program main
   USE base_var
   IMPLICIT NONE
   integer::ip,interval,nend,nstart,ntime,istep,nstep
   character(len=3)::ctmp3
   character(len=1024)::ctmp
   character(len=2048)::cline

 100 continue
   write(*,*)"Input number of processor"
   read(*,*,err=100,end=100)np
   if(np==0)then
     write(*,*)"Number of processor should be larger than 0"
     goto 100
   endif

 106 continue
   write(*,*)"Input number of digit for time step"
   read(*,*,err=106,end=106)ndigit
   if(ndigit==0)then
     write(*,*)"Number of digit should be larger than 0"
     goto 106
   endif

 101 continue
   write(*,*)"Input start of time step"
   read(*,*,err=101,end=101)nstart

 102 continue
   write(*,*)"Input end of time step"
   read(*,*,err=102,end=102)nend
   if(nstart>nend)then
     write(*,*)"nstart should be smaller than nend"
     goto 101
   endif

 103 continue
   write(*,*)"Input time step increment"
   read(*,*,err=103,end=103)interval
 104 continue
   write(*,*)"Input common file name."
   write(*,*)"Example: uvw-c-press_p"
   read(*,*,err=104,end=104)fncommon
 105 continue
   write(*,*)"Input number of variables."
   read(*,*,err=105,end=105)nvariable


   nstep=1+(nend-nstart)/interval

   allocate(fname(np))

   DO istep=1,nstep
      ntime=nstart+(istep-1)*interval

      call int2char(ntime,ctmp,ndigit)
      DO ip=1,np
         call int2char3(ip-1,ctmp3)
         fname(ip)=trim(fncommon)//ctmp3//"_"//ctmp(1:ndigit)//".dat"
         write(*,*)trim(fname(ip))
      ENDDO

      IF(istep==1)call alloc
      call gather
      call output(ntime)
      call delfile
   ENDDO

   deallocate(x,v)

   stop
   end
!
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
!
!-----------------------------------------------------------------------
   SUBROUTINE INT2CHAR3( NUM, CWRK )
   INTEGER NUM, J
   CHARACTER(LEN=3) CWRK
!
   CWRK = '   '
   WRITE( CWRK, '( I3 )' ) NUM
   DO J = 1, 3
       IF ( CWRK( J:J ).EQ.' ' ) CWRK( J:J ) = '0'
   END DO
   RETURN
   END SUBROUTINE INT2CHAR3
!
!-----------------------------------------------------------------------
   SUBROUTINE alloc
   USE base_var
   IMPLICIT NONE
   INTEGER::i,j,k,ifile,idummy
   CHARACTER(len=12)::cdummy1,cdummy2

   ! set icmax,jcmax,kcmax --------------------------
   icmax=0
   jcmax=0
   kcmax=0
   DO ifile=1,np
      OPEN(10,file=fname(ifile),status='OLD',err=999)
         DO i=1,100
            READ(10,*,err=998,end=998)cdummy1,cdummy2
            IF(cdummy2=="I-RANGE")EXIT
         ENDDO
         BACKSPACE(10)
         READ(10,*,err=998,end=998)cdummy1,cdummy2,idummy,i
         READ(10,*,err=998,end=998)cdummy1,cdummy2,idummy,j
         READ(10,*,err=998,end=998)cdummy1,cdummy2,idummy,k
         IF(icmax<=i)icmax=i
         IF(jcmax<=j)jcmax=j
         IF(kcmax<=k)kcmax=k
      CLOSE(10)
   ENDDO
   WRITE(*,*)"icmax,jcmax,kcmax=",icmax,jcmax,kcmax
   inmax=icmax+1
   jnmax=jcmax+1
   knmax=kcmax+1

   ALLOCATE(x(inmax,jnmax,knmax,3))
   ALLOCATE(v(icmax,jcmax,kcmax,nvariable))

   RETURN
 999 continue
   WRITE(*,*)"Error!!! File not found. ",trim(fname(ifile))
   STOP
 998 continue
   WRITE(*,*)"Error!!! File read error. ",trim(fname(ifile))
   STOP

   END SUBROUTINE alloc

!-----------------------------------------------------------------------
   SUBROUTINE gather
   USE base_var
   IMPLICIT NONE
   INTEGER::ics,ice,jcs,jce,kcs,kce   !!cell
   INTEGER::ins,ine,jns,jne,kns,kne   !!node
   INTEGER::i,j,k,ifile,iline,idummy,m
   CHARACTER(len=12)::cdummy1,cdummy2

   DO ifile=1,np
      OPEN(10,file=fname(ifile),status='OLD',err=999)
         DO i=1,100
            READ(10,*,err=998,end=998)cdummy1,cdummy2
            IF(cdummy2=="I-RANGE")EXIT
         ENDDO
         BACKSPACE(10)
         READ(10,*,err=998,end=998)cdummy1,cdummy2,ics,ice
         READ(10,*,err=998,end=998)cdummy1,cdummy2,jcs,jce
         READ(10,*,err=998,end=998)cdummy1,cdummy2,kcs,kce
         ins=ics; ine=ice+1
         jns=jcs; jne=jce+1
         kns=kcs; kne=kce+1
         READ(10,*,err=998,end=998) !!VARIABLES
         READ(10,*,err=998,end=998) !!ZONE
         DO m=1,3
            READ(10,*,err=998,end=998) !!# COORDINATES
            READ(10,*,err=998,end=998)  &
              (((x(i,j,k,m),i=ins,ine),j=jns,jne),k=kns,kne)
         ENDDO
         DO m=1,nvariable
            READ(10,*,err=998,end=998) !!# VARIABLES
            READ(10,*,err=998,end=998)  &
              (((v(i,j,k,m),i=ics,ice),j=jcs,jce),k=kcs,kce)
         ENDDO
      CLOSE(10)
   ENDDO

   RETURN
 999 continue
   WRITE(*,*)"Error!!! File not found. ",trim(fname(ifile))
   STOP
 998 continue
   WRITE(*,*)"Error!!! File read error. ",trim(fname(ifile))
   STOP

   END SUBROUTINE gather
!-----------------------------------------------------------------------
   SUBROUTINE output(nt)
   USE base_var
   IMPLICIT NONE
   INTEGER::nt
   INTEGER::i,j,k,ifile,iline,idummy,m
   CHARACTER(len=2048)::fout,cline
   CHARACTER(len=12)::cdummy1,cdummy2
   CHARACTER(len=1024)::ctmp

   call int2char(nt,ctmp,ndigit)
   fout=trim(fncommon)//"all_"//ctmp(1:ndigit)//".dat"
   write(*,*)"Output to ",trim(fout)

   OPEN(10,file=fname(1),status='OLD')
   OPEN(11,file=fout,status='UNKNOWN')
      !DO iline=1,4
      !   READ(10,*)
      !ENDDO
      DO i=1,100
         READ(10,*)cdummy1,cdummy2
         IF(cdummy2=="K-RANGE")EXIT
      ENDDO

      READ(10,'(a2048)')cline  !!VARIABLES
      WRITE(11,*)trim(cline)
      WRITE(11,'(a,i9,a,i9,a,i9)')" ZONE I=",inmax," J=",jnmax," k=",knmax
      IF(nvariable==1)THEN
         WRITE(11,'(a)')"DATAPACKING=BLOCK, VARLOCATION=([4]=CELLCENTERED)"
      ELSE IF(nvariable+3<=9)THEN
         WRITE(11,'(a,i1,a)')"DATAPACKING=BLOCK, VARLOCATION=([4-"  &
                 ,nvariable+3,"]=CELLCENTERED)"
      ELSE
         WRITE(11,'(a,i2,a)')"DATAPACKING=BLOCK, VARLOCATION=([4-"  &
                 ,nvariable+3,"]=CELLCENTERED)"
      ENDIF
      DO m=1,3
         WRITE(11,'(a,i)')"# COORDINATE-",m
         WRITE(11,*)  &
           (((x(i,j,k,m),i=1,inmax),j=1,jnmax),k=1,knmax)
      ENDDO
      DO m=1,nvariable
         WRITE(11,'(a10,i1)')"# VARIABALE-",m
         WRITE(11,*)  &
           (((v(i,j,k,m),i=1,icmax),j=1,jcmax),k=1,kcmax)
      ENDDO
   CLOSE(10)
   CLOSE(11)

   RETURN

   END SUBROUTINE output

!-----------------------------------------------------------------------
   SUBROUTINE delfile
   USE base_var
   IMPLICIT NONE
   INTEGER::ifile

   DO ifile=1,np
      OPEN(10,file=fname(ifile),status='OLD')
      CLOSE(10,status='delete')
   ENDDO

   RETURN
   END SUBROUTINE delfile

