!  generate mcr(script) file for Tecplot
   IMPLICIT NONE
   integer::ip,interval,nend,ntime,istep,nstep,ndigit,nval,i,itmp
   character(len=1024)::ctmp,vlist
   character(len=2048)::fncommon,cline,fname,fout

 101 continue
   write(*,*)"Input start of time step"
   read(*,*,err=101,end=101)ntime
 102 continue
   write(*,*)"Input end of time step"
   read(*,*,err=102,end=102)nend
 103 continue
   write(*,*)"Input time step increment"
   read(*,*,err=103,end=103)interval
 104 continue
   write(*,*)"Input common file name."
   write(*,*)"Example: uvw-c-press_pall_"
   read(*,*,err=104,end=104)fncommon
 106 continue
   write(*,*)"Input number of digit for time step"
   read(*,*,err=106,end=106)ndigit
   if(ndigit==0)then
     write(*,*)"Number of digit should be larger than 0"
     goto 106
   endif
 108 continue
   write(*,*)"Input number of variables. For example, uvw-c is 4"
   read(*,*,err=108,end=108)nval
   if(nval<0) goto 108
 105 continue
   write(*,*)"Input output-file name.  Example:animation.mcr"
   read(*,*,err=105,end=105)fout

!--make vlist
   vlist="X Y Z U V W"
   do i=1,nval-3
     !write(*,*)char(64+i)
     vlist=trim(vlist)//" "//char(64+i)
     !write(*,*)trim(vlist)
   enddo

!--output 

   OPEN(10,file=fout,status="unknown")

   nstep=(nend-ntime)/interval+1

   call int2char(ntime,ctmp,ndigit)
   fname=trim(fncommon)//ctmp(1:ndigit)//".plt"

   write(10,'(a)')"#!MC 1100"
   cline="$!READDATASET  '"
   cline=trim(cline)//" "//trim(fname)
   cline=trim(cline)//"'"
   write(10,'(a)')trim(cline)
   write(10,'(a)')"  READDATAOPTION = NEW"
   write(10,'(a)')"  RESETSTYLE = NO"
   write(10,'(a)')"  INCLUDETEXT = NO"
   write(10,'(a)')"  INCLUDEGEOM = NO"
   write(10,'(a)')"  INCLUDECUSTOMLABELS = NO"
   write(10,'(a)')"  VARLOADMODE = BYNAME"
   write(10,'(a)')"  ASSIGNSTRANDIDS = YES"
   write(10,'(a)')"  INITIALPLOTTYPE = CARTESIAN3D"
   write(10,'(a)')"$!EXPORTSETUP EXPORTFORMAT = PNG"
   write(10,'(a)')"$!EXPORTSETUP IMAGEWIDTH = 1200"
   write(10,'(3a)')"$!EXPORTSETUP EXPORTFNAME = './",ctmp(1:ndigit),".png'"
   write(10,'(a)')"$!EXPORT"
   write(10,'(a)')"  EXPORTREGION = CURRENTFRAME"
   write(10,*)

   DO istep=1,nstep-1
      ntime=ntime+interval

      call int2char(ntime,ctmp,ndigit)
      fname=trim(fncommon)//ctmp(1:ndigit)//".plt"

      cline="$!READDATASET  '"
      cline=trim(cline)//" "//trim(fname)
      cline=trim(cline)//"'"
      write(10,'(a)')trim(cline)

      write(10,*)"READDATAOPTION = NEW"
      write(10,*)"  RESETSTYLE = NO"
      write(10,*)"  INCLUDETEXT = NO"
      write(10,*)"  INCLUDEGEOM = NO"
      write(10,*)"  INCLUDECUSTOMLABELS = NO"
      write(10,*)"  VARLOADMODE = BYNAME"
      write(10,*)"  ASSIGNSTRANDIDS = YES"
      write(10,'(a,a,a)')"  VARNAMELIST = ' ",trim(vlist)," '"
      write(10,*)"$!REDRAW "
      write(10,'(3a)')"$!EXPORTSETUP EXPORTFNAME = './",ctmp(1:ndigit),".png'"
      write(10,'(a)')"$!EXPORT"
      write(10,'(a)')"  EXPORTREGION = CURRENTFRAME"
      write(10,*)
   ENDDO

   write(10,*)"$!EXPORTFINISH "

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
/*-----------------------------------------------------------------------------+
 '$Id: tecmcr-png.f90,v 1.2 2018/09/17 12:53:28 sato Exp $'/
+-----------------------------------------------------------------------------*/
