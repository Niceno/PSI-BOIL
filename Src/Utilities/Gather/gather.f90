!  gather tecplot (PSI-Boil dat format) file
   module base_var
     IMPLICIT NONE
     REAL(4),ALLOCATABLE::x(:,:,:,:),v(:,:,:,:)
     INTEGER::icmax,jcmax,kcmax   !!cell
     INTEGER::inmax,jnmax,knmax   !!node
     INTEGER::np,nvariable,ndigit
     INTEGER::nodal
     LOGICAL::withBuffer
     CHARACTER(len=2048),allocatable::fname(:)
     CHARACTER(len=24),allocatable::valname(:)
     CHARACTER(len=2048)fncommon
     CHARACTER(len=2048)fname_out
   end module base_var

!-----------------------------------------------------------------------
 program main
   USE base_var
   IMPLICIT NONE
   integer::ip,interval,nend,nstart,ntime,istep,nstep
   character(len=3)::ctmp3
   character(len=4)::ctmp4
   character(len=5)::ctmp5
   character(len=1024)::ctmp
   character(len=2048)::cline

 100 continue
   write(*,*)"Input number of processor"
   read(*,*,err=100,end=100)np
   if(np==0)then
     write(*,*)"Number of processor should be larger than 0"
     !goto 100
     stop
   endif

 106 continue
   write(*,*)"Input number of digit for time step"
   read(*,*,err=106,end=106)ndigit
   if(ndigit==0)then
     write(*,*)"Number of digit should be larger than 0"
     !goto 106
     stop
   endif

 101 continue
   write(*,*)"Input start of time step"
   read(*,*,err=101,end=101)nstart

 102 continue
   write(*,*)"Input end of time step"
   read(*,*,err=102,end=102)nend
   if(nstart>nend)then
     write(*,*)"nstart should be smaller than nend"
     !goto 101
     stop
   endif

 103 continue
   write(*,*)"Input time step increment"
   read(*,*,err=103,end=103)interval
 104 continue
   write(*,*)"Input common file name."
   write(*,*)"Example: uvw-c-press_p"
   read(*,*,err=104,end=104)fncommon
! 105 continue
!   write(*,*)"Input number of variables."
!   read(*,*,err=105,end=105)nvariable


   nstep=1+(nend-nstart)/interval

   allocate(fname(np))

   DO istep=1,nstep
      ntime=nstart+(istep-1)*interval

      call int2char(ntime,ctmp,ndigit)
      DO ip=1,np
         if(ip-1<=999)then
           call int2char3(ip-1,ctmp3)
           fname(ip)=trim(fncommon)//ctmp3//"_"//ctmp(1:ndigit)//".dat"
         else if(ip-1<=9999)then
           call int2char4(ip-1,ctmp4)
           fname(ip)=trim(fncommon)//ctmp4//"_"//ctmp(1:ndigit)//".dat"
         else
           call int2char5(ip-1,ctmp5)
           fname(ip)=trim(fncommon)//ctmp5//"_"//ctmp(1:ndigit)//".dat"
         endif
         write(*,*)trim(fname(ip))
      ENDDO
      fname_out=trim(fncommon)//"all_"//ctmp(1:ndigit)//".plt"
#ifdef DEBUG
      WRITE(*,*)'call getnvar'
#endif
      call getnvar

#ifdef DEBUG
      WRITE(*,*)'call alloc'
#endif
      call alloc

#ifdef DEBUG
      WRITE(*,*)'call gather'
#endif
      call gather
      call output(ntime)
#ifdef ZIP
      call compress
#endif
#ifndef SZPLT
      call delfile
#endif
      call dealloc
   ENDDO

   stop
   end
!
!-----------------------------------------------------------------------
#ifdef ZIP
   SUBROUTINE compress
   USE base_var
   CHARACTER(len=2048)::command
   command="gzip "//trim(fname_out)
      write(*,*)'command= ',trim(command)
      call system(trim(command))

   RETURN

   END SUBROUTINE compress
#endif
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
   SUBROUTINE INT2CHAR4( NUM, CWRK )
   INTEGER NUM, J
   CHARACTER(LEN=4) CWRK
!
   CWRK = '   '
   WRITE( CWRK, '( I4 )' ) NUM
   DO J = 1, 4
       IF ( CWRK( J:J ).EQ.' ' ) CWRK( J:J ) = '0'
   END DO
   RETURN
   END SUBROUTINE INT2CHAR4
!
!-----------------------------------------------------------------------
   SUBROUTINE INT2CHAR5( NUM, CWRK )
   INTEGER NUM, J
   CHARACTER(LEN=5) CWRK
!
   CWRK = '   '
   WRITE( CWRK, '( I5 )' ) NUM
   DO J = 1, 5
       IF ( CWRK( J:J ).EQ.' ' ) CWRK( J:J ) = '0'
   END DO
   RETURN
   END SUBROUTINE INT2CHAR5
!
!-----------------------------------------------------------------------
   SUBROUTINE getnvar
   USE base_var
   IMPLICIT NONE
   INTEGER::i,j,k,ifile,idummy,m
   CHARACTER(len=12)::cdummy1,cdummy2
   CHARACTER(len=2048)::cline1,cline2
#ifdef DEBUG
   WRITE(*,*)'getnvar:start'
   WRITE(*,*)'getnvar:fname(1)= ',trim(fname(1))
#endif

   ! set nvariable ----------------------------------
   OPEN(10,file=trim(fname(1)),status='OLD',err=999)
   !OPEN(10,file='uvw-c-tpr-press-mdot_p000_000002.dat',status='OLD',err=999)
#ifdef DEBUG
      WRITE(*,*)'getnvar:open file= ',trim(fname(1))
#endif
      DO i=1,100
         READ(10,'(a2048)')cline1
         IF(cline1(1:9)=="VARIABLES" .OR. cline1(2:10)=="VARIABLES") EXIT
      ENDDO
   CLOSE(10)

   cline2=""
   cline2=cline1(11:2048)
#ifdef DEBUG
   WRITE(*,*)'getnvar:cline2= ',trim(cline2)
#endif

   j=0
   DO i=1,2048
      if(cline2(i:i)=='"') j=j+1
   ENDDO
   nvariable=j/2-3
   ALLOCATE(valname(nvariable+3))
   valname=""

   j=0
   k=0
   m=0
   DO i=1,2048
      IF (cline2(i:i)=='"'.AND.m==0) THEN
         j=j+1
         m=1
         k=0
      ELSE IF(cline2(i:i)=='"'.AND.m==1) THEN
         m=0
         IF(j==nvariable+3) EXIT
      ELSE IF(j>=1) THEN
        k=k+1
        valname(j)(k:k)=cline2(i:i)
      ENDIF
   ENDDO

   DO m=1,nvariable+3
      WRITE(*,*)m,trim(valname(m))
   ENDDO

   RETURN

 999 continue
   WRITE(*,*)"getnvar:Error!!! File not found. ",trim(fname(ifile))
   STOP
 998 continue
   WRITE(*,*)"getnvar:Error!!! File read error. ",trim(fname(ifile))
   STOP

   END

!-----------------------------------------------------------------------
   SUBROUTINE alloc
   USE base_var
   IMPLICIT NONE
   INTEGER::i,j,k,ifile,idummy,m
   CHARACTER(len=12)::cdummy1,cdummy2,cdummy3,cdummy4,cdummy5

   ! set icmax,jcmax,kcmax --------------------------
   icmax=0
   jcmax=0
   kcmax=0
   withBuffer=.FALSE.
   DO ifile=1,np
      OPEN(11,file=fname(ifile),status='OLD',err=999)
         READ(11,*,err=998,end=998)
         READ(11,*,err=500,end=500)cdummy1,cdummy2,cdummy3,cdummy4,cdummy5
         IF(trim(cdummy5)=="BUFFER")THEN
           withBuffer=.TRUE.
           goto 501
         ENDIF

 500     CONTINUE
         BACKSPACE(11)
         BACKSPACE(11)
         READ(11,*,err=998,end=998)cdummy1,cdummy2
         IF(cdummy2=="NODAL")THEN
           nodal=1
         ELSE
           nodal=0
         ENDIF

 501     CONTINUE

         READ(11,*,err=998,end=998)cdummy1,cdummy2,idummy,i
         READ(11,*,err=998,end=998)cdummy1,cdummy2,idummy,j
         READ(11,*,err=998,end=998)cdummy1,cdummy2,idummy,k
         IF(icmax<=i)icmax=i
         IF(jcmax<=j)jcmax=j
         IF(kcmax<=k)kcmax=k
      CLOSE(11)
   ENDDO
   WRITE(*,*)'WITH BUFFER= ',withBuffer
   IF(withBuffer)THEN
     icmax=icmax+2
     jcmax=jcmax+2
     kcmax=kcmax+2
   ENDIF
   WRITE(*,*)"icmax,jcmax,kcmax=",icmax,jcmax,kcmax
   inmax=icmax+1
   jnmax=jcmax+1
   knmax=kcmax+1

   ALLOCATE(x(inmax,jnmax,knmax,3))
   IF (nodal==0) THEN
     ALLOCATE(v(icmax,jcmax,kcmax,nvariable))
   ELSE
     ALLOCATE(v(inmax,jnmax,knmax,nvariable))
   ENDIF

   RETURN
 999 continue
   WRITE(*,*)"alloc:Error!!! File not found. ",trim(fname(ifile))
   STOP
 998 continue
   WRITE(*,*)"alloc:Error!!! File read error. ",trim(fname(ifile))
   STOP

   END SUBROUTINE alloc

!-----------------------------------------------------------------------
   SUBROUTINE dealloc
   USE base_var
   IMPLICIT NONE
     DEALLOCATE(x,v)
     DEALLOCATE(valname)
   RETURN
   END SUBROUTINE dealloc
!-----------------------------------------------------------------------
   SUBROUTINE gather
   USE base_var
   IMPLICIT NONE
   INTEGER::ics,ice,jcs,jce,kcs,kce   !!cell
   INTEGER::ins,ine,jns,jne,kns,kne   !!node
   INTEGER::i,j,k,ifile,iline,idummy,m
   CHARACTER(len=12)::cdummy1,cdummy2

   x=1.0d+20

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

         IF(withBuffer)THEN
           ins=ics
           ine=ice+3
           jns=jcs
           jne=jce+3
           kns=kcs
           kne=kce+3

           ice=ice+2
           jce=jce+2
           kce=kce+2
         ELSE
           ins=ics; ine=ice+1
           jns=jcs; jne=jce+1
           kns=kcs; kne=kce+1
         ENDIF

         !WRITE(*,*)ins,ine,jns,jne,kns,kne
         READ(10,*,err=998,end=998) !!VARIABLES
         READ(10,*,err=998,end=998) !!ZONE
         DO m=1,3
            !WRITE(*,*)'m=',m
            READ(10,*,err=998,end=998) !!# COORDINATES
            READ(10,*,err=998,end=998)  &
              (((x(i,j,k,m),i=ins,ine),j=jns,jne),k=kns,kne)
         ENDDO
         IF(nodal==0)THEN
           DO m=1,nvariable
              READ(10,*,err=998,end=998) !!# VARIABLES
              READ(10,*,err=998,end=998)  &
                (((v(i,j,k,m),i=ics,ice),j=jcs,jce),k=kcs,kce)
           ENDDO
         ELSE
           DO m=1,nvariable
              READ(10,*,err=998,end=998) !!# VARIABLES
              READ(10,*,err=998,end=998)  &
                (((v(i,j,k,m),i=ins,ine),j=jns,jne),k=kns,kne)
           ENDDO
         ENDIF
      CLOSE(10)
   ENDDO

   !Check
   IF (maxval(x)>1.0d+20) THEN
     WRITE(*,*)"gather:Error!!! Did not read all files."
     STOP
   ENDIF

   RETURN
 999 continue
   WRITE(*,*)"gather:Error!!! File not found. ",trim(fname(ifile))
   STOP
 998 continue
   WRITE(*,*)"gather:Error!!! File read error. ",trim(fname(ifile))
   STOP

   END SUBROUTINE gather
!-----------------------------------------------------------------------
   SUBROUTINE output(nt)
   USE base_var
   IMPLICIT NONE
   INCLUDE 'tecio.f90'
   INTEGER::nt
   INTEGER::i,j,k,ifile,iline,idummy,m
   CHARACTER(len=128)::fout
   CHARACTER(len=2048)::cline
   CHARACTER(len=12)::cdummy1,cdummy2
   CHARACTER(len=1024)::ctmp
   REAL(4),ALLOCATABLE::anode(:,:,:),acell(:,:,:)

   character*1 NULLCHR
   Integer*4   Debug,III,IIII,NPts,NElm
   Real*8    SolTime
   Integer*4 VIsDouble, FileType, FileFormat
   Integer*4 ZoneType,StrandID,unused,IsBlock
   Integer*4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
   Integer*4 Valuelocation(nvariable+3)
   POINTER   (NullPtr,Null)
   Integer*4 Null(*)

   NULLCHR = CHAR(0)
   NullPtr = 0
   Debug   = 0
   FileType = 0
#ifndef SZPLT
   FileFormat = 0 ! 0 = PLT, 1 = SZPLT
#else
   FileFormat = 1 ! 0 = PLT, 1 = SZPLT
#endif
   VIsDouble = 0
   ZoneType = 0
   SolTime = 0.0
   StrandID = 0
   unused = 0 ! ParentZone is no longer used
   IsBlock = 1
   ICellMax = 0
   JCellMax = 0
   KCellMax = 0
   NFConns = 0
   FNMode = 0
   ShrConn = 0

   DO i=1,3
     Valuelocation(i)=1
   ENDDO
   IF (nodal==0) THEN
     DO i=4,3+nvariable
       Valuelocation(i)=0
     ENDDO
   ELSE
     DO i=4,3+nvariable
       Valuelocation(i)=1
     ENDDO
   ENDIF

   ALLOCATE(anode(inmax,jnmax,knmax),acell(inmax-1,jnmax-1,knmax-1))
!
!... Set output file name
!
   call int2char(nt,ctmp,ndigit)
#ifndef SZPLT
   fout=trim(fncommon)//"all_"//ctmp(1:ndigit)//".plt"
#else
   fout=trim(fncommon)//"all_"//ctmp(1:ndigit)
#endif
   write(*,*)"Output to ",trim(fout)
!
!... Set variable name
!
   cline=valname(1)
   DO m=2,nvariable+3
     cline=trim(cline)//" "//valname(m)
   ENDDO
   !WRITE(*,*)trim(cline)
!
!... Open the file and write the tecplot datafile 
!... header information.
!
   I = TecIni142('DATASET'//NULLCHR, &
                 trim(cline)//NULLCHR, &
                 trim(fout)//NULLCHR, &
                 '.'//NULLCHR, &
                 FileFormat, &
                 FileType, &
                 Debug, &
                 VIsDouble)
!
!... Write the zone header information.
!
#ifndef VISIT
   I = TecZne142(ctmp(1:ndigit)//NULLCHR, &
#else
   I = TecZne142('000000'//NULLCHR, &
#endif
                 ZoneType, &
                 inmax, &
                 jnmax, &
                 knmax, &
                 ICellMax, &
                 JCellMax, &
                 KCellMax, &
                 SolTime, &
                 StrandID, &
                 unused, &
                 IsBlock, &
                 NFConns, &
                 FNMode, &
                 0, &
                 0, &
                 0, &
                 Null, &
                 Valuelocation, &
                 Null, &
                 ShrConn)
!
!... Write out the field data.
!
   III = inmax*jnmax*knmax
   IIII = (inmax-1)*(jnmax-1)*(knmax-1)
   DO m=1,3
     DO k=1,knmax
     DO j=1,jnmax
     DO i=1,inmax
       anode(i,j,k)=x(i,j,k,m)
     ENDDO
     ENDDO
     ENDDO
     I   = TecDat142(III,anode,0)
   ENDDO

   IF(nodal==0)THEN
     DO m=1,nvariable
       DO k=1,knmax-1
       DO j=1,jnmax-1
       DO i=1,inmax-1
         acell(i,j,k)=v(i,j,k,m)
       ENDDO
       ENDDO
       ENDDO
       I = TecDat142(IIII,acell,0)
     ENDDO
   ELSE
     DO m=1,nvariable
       DO k=1,knmax
       DO j=1,jnmax
       DO i=1,inmax
         anode(i,j,k)=v(i,j,k,m)
       ENDDO
       ENDDO
       ENDDO
       I = TecDat142(III,anode,0)
     ENDDO
   ENDIF

   I = TecEnd142()

   DEALLOCATE(anode,acell)

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

