!  generate mcr(script) file for Tecplot
   IMPLICIT NONE
   integer::ip,interval,nend,ntime,istep,nstep,ndigit,ispeed,nval,i,itmp,i_png_avi,i_width
   integer::nplt_each
   character(len=1024)::ctmp,vlist
   character(len=2048)::fncommon,fncommon2,cline,fname,fout,fname2

 100 continue
   write(*,*)"Input number of plt-files drawn in each picture. 1, 2 or 3"
   read(*,*,err=100,end=100)nplt_each
   if (nplt_each/=1.and.nplt_each/=2.and.nplt_each/=3) goto 100
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
!106 continue
!  write(*,*)"Input number of digit for time step"
!  read(*,*,err=106,end=106)ndigit
!  if(ndigit==0)then
!    write(*,*)"Number of digit should be larger than 0"
!    goto 106
!  endif
   ndigit=6
!107 continue
!  write(*,*)"Input number of variables. For example, uvw-c is 4"
!  read(*,*,err=107,end=107)nval
!  if(nval<0) goto 107
 108 continue
   if(nplt_each==2)then
     write(*,*)"Input the next file name. The name must be fixed. Ex: nozzle.plt"
     read(*,*,err=104,end=104)fname2
   endif
   if(nplt_each==3)then
     write(*,*)"Input common file name of the next file. Example: particles_"
     read(*,*,err=104,end=104)fncommon2

     write(*,*)"Input the next file name. The name must be fixed. Ex: nozzle.plt"
     read(*,*,err=104,end=104)fname2
   endif

!--output setting
 110 continue
   write(*,*)"Input output-file type. 1:PNG, 2:AVI"
   read(*,*,err=110,end=110)i_png_avi
   if (i_png_avi/=1.and.i_png_avi/=2) goto 110
 111 continue
   write(*,*)"Input output-file name.  Example:png.mcr, avi.mcr"
   read(*,*,err=111,end=111)fout
 112 continue
   write(*,*)"Input output picture width.  Example:1280"
   read(*,*,err=112,end=112)i_width

   if(i_png_avi==2)then
 113 continue
     write(*,*)"Input animation speed (frame/sec.)"
     read(*,*,err=113,end=113)ispeed
   endif

!--make vlist
!  vlist="X Y Z U V W"
!  do i=1,nval-3
!    vlist=trim(vlist)//" "//char(64+i)
!  enddo

!--output 

   OPEN(10,file=fout,status="unknown")

   nstep=(nend-ntime)/interval+1

   call int2char(ntime,ctmp,ndigit)
   fname=trim(fncommon)//ctmp(1:ndigit)//".plt"
   WRITE(*,*)trim(fname)

   write(10,'(a)')"#!MC 1100"
   cline="$!READDATASET  '"
   cline=trim(cline)//trim(fname)
   if(nplt_each==2)then
     cline=trim(cline)//" "//trim(fname2)
   endif
   if(nplt_each==3)then
     WRITE(*,*)trim(cline)
     fname=trim(fncommon2)//ctmp(1:ndigit)//".plt"
     WRITE(*,*)trim(fname)
     cline=trim(cline)//" "//trim(fname)
     WRITE(*,*)trim(cline)
     cline=trim(cline)//" "//trim(fname2)
     WRITE(*,*)trim(cline)
   endif

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
   if(i_png_avi==1)then
     write(10,'(a)')"$!EXPORTSETUP EXPORTFORMAT = PNG"
     write(10,'(a,i4)')"$!EXPORTSETUP IMAGEWIDTH = ",i_width
     write(10,'(3a)')"$!EXPORTSETUP EXPORTFNAME = './",ctmp(1:ndigit),".png'"
     write(10,'(a)')"$!EXPORT"
   else if(i_png_avi==2)then
     write(10,'(a)')"$!EXPORTSETUP EXPORTFORMAT = AVI"
     write(10,'(a,i4)')"$!EXPORTSETUP IMAGEWIDTH = ",i_width
     write(10,'(a,i4)')"$!EXPORTSETUP ANIMATIONSPEED = ", ispeed
     write(10,'(a)')"$!EXPORTSETUP EXPORTFNAME = './export.avi'"
     write(10,'(a)')"$!EXPORTSTART "
     write(10,'(a)')"  EXPORTREGION = CURRENTFRAME"
   endif
   write(10,*)

   DO istep=1,nstep-1
      ntime=ntime+interval

      call int2char(ntime,ctmp,ndigit)
      fname=trim(fncommon)//ctmp(1:ndigit)//".plt"

      cline="$!READDATASET  '"
      cline=trim(cline)//trim(fname)
      if(nplt_each==2)then
        cline=trim(cline)//" "//trim(fname2)
      endif
      if(nplt_each==3)then
        WRITE(*,*)trim(cline)
        fname=trim(fncommon2)//ctmp(1:ndigit)//".plt"
        WRITE(*,*)trim(fname)
        cline=trim(cline)//" "//trim(fname)
        WRITE(*,*)trim(cline)
        cline=trim(cline)//" "//trim(fname2)
        WRITE(*,*)trim(cline)
      endif
      cline=trim(cline)//"'"
      write(10,'(a)')trim(cline)

      write(10,*)"READDATAOPTION = NEW"
      write(10,*)"  RESETSTYLE = NO"
      write(10,*)"  INCLUDETEXT = NO"
      write(10,*)"  INCLUDEGEOM = NO"
      write(10,*)"  INCLUDECUSTOMLABELS = NO"
      write(10,*)"  VARLOADMODE = BYNAME"
      write(10,*)"  ASSIGNSTRANDIDS = YES"
      !write(10,*)"  VARNAMELIST = 'X Y Z U V W A B'"
      !write(10,'(a,a,a)')"  VARNAMELIST = ' ",trim(vlist)," '"
      write(10,*)"$!REDRAW "
      if(i_png_avi==1)then
        write(10,'(3a)')"$!EXPORTSETUP EXPORTFNAME = './",ctmp(1:ndigit),".png'"
        write(10,'(a)')"$!EXPORT"
        write(10,'(a)')"  EXPORTREGION = CURRENTFRAME"
      else if(i_png_avi==2)then
        write(10,*)"$!EXPORTNEXTFRAME "
      endif
      write(10,*)
   ENDDO

   if(i_png_avi==2)then
     write(10,*)"$!EXPORTFINISH "
   endif

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
