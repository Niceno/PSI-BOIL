  Program pathline
  implicit none
  integer::nstart,nend,interval,nstep,istep,itime,ndigit,nval,np,n_header
  integer::iline,i,j,num_words
  real(8),allocatable::x(:,:),y(:,:),z(:,:),u(:,:),v(:,:),w(:,:),time(:)
  real(8),allocatable::val(:,:,:)  ! ival,ip,istep
  character(len=512),allocatable::fname(:),vname(:)
  character(len=512)::fncommon
  character(len=512)::ctmp,ctmp1
  character(len=20 )::format_str
  character(len=200),allocatable::words(:)

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

  ! set np, nval
  nval = 0
  OPEN(10,file=trim(fname(1)),status='old',err=999)

    ! nval
    DO iline = 1,10
      READ(10,*)ctmp
      if (trim(ctmp) .eq. "VARIABLES=") then
        BACKSPACE(10)
        READ(10,'(512a)')ctmp
        WRITE(*,*)trim(ctmp)
        CALL split_string(ctmp, num_words)
        nval = num_words-1
        EXIT
      endif
    ENDDO
    WRITE(*,*)'nval=',nval
    IF (nval==0) THEN
      WRITE(*,*)'ERROR! nval=0'
      STOP
    ENDIF
    ALLOCATE(vname(nval))

    ! np
    REWIND(10)
    DO iline = 1,10
      READ(10,*)ctmp
      !WRITE(*,*)trim(ctmp)
      if (trim(ctmp) .eq. "ZONE") then
        n_header = iline
        BACKSPACE(10)
        READ(10,*)ctmp,ctmp,np
        EXIT
      endif
    ENDDO
    WRITE(*,*)'np=',np

    ! vname
    REWIND(10)
    DO iline = 1,10
      READ(10,*)ctmp
      if (trim(ctmp) .eq. "VARIABLES=") then
        BACKSPACE(10)
        READ(10,*)ctmp,(vname(j),j=1,nval)
        EXIT
      endif
    ENDDO
    WRITE(*,*)'VARIABLES=',(' ',trim(vname(j)),j=1,nval)
  CLOSE(10)

  ! allocate variables
  allocate(x(np,nstep),y(np,nstep),z(np,nstep), &
           u(np,nstep),v(np,nstep),w(np,nstep),time(nstep))
  allocate(val(nval-6,np,nstep))

  ! read data
  DO istep = 1, nstep
    OPEN(10,file=trim(fname(istep)),status='old',err=998)
    !WRITE(*,*)'n_header=',n_header
    DO iline = 1, n_header
      READ(10,*)
    ENDDO
    DO i = 1, np
      READ(10,*)x(i,istep),y(i,istep),z(i,istep), &
                u(i,istep),v(i,istep),w(i,istep), &
                (val(j,i,istep),j=1,nval-6)
    ENDDO
    CLOSE(10)
    WRITE(*,*)trim(fname(istep))
  ENDDO

  OPEN(10,file='pathlines.dat',FORM='formatted')
    WRITE(*,*)'# output to: pathlines.dat'
    WRITE(format_str, '(A,I0,A)') '(', nval, 'E16.8)'
    DO i = 1, np
      call int2char(i,ctmp1,ndigit)
      ctmp = 'TITLE= "Particle'//ctmp1(1:ndigit)//'"'
      !WRITE(10,*)'TITLE= "Particle',i,'"'
      WRITE(10,*)trim(ctmp)
      WRITE(10,*)'VARIABLES=',(' ',trim(vname(j)),j=1,nval)
      !WRITE(10,*)"VARIABLES= X Y Z U V W"
      WRITE(10,*)"ZONE I=",nstep, "DATAPACKING=POINT"

      !DO i = 1, nval
      !  WRITE(format_str, '(A,I0,A)') '(', i, 'E16.8)'
      !  write(*, format_str) 1.0  ! ここで1.0は仮の値です
      !ENDDO
      !WRITE(format_str, '(A,I0,A)') '(', nval, 'E16.8)'
      !WRITE(*,*)'format_str=',format_str
      DO istep = 1, nstep
        WRITE(10,format_str)x(i,istep),y(i,istep),z(i,istep), &
                u(i,istep),v(i,istep),w(i,istep),  &
                (val(j,i,istep),j=1,nval-6)
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

!-----------------------------------------------------------------------
   SUBROUTINE split_string(input_line, num_words)
        character(len=*), intent(in) :: input_line
        integer, intent(out) :: num_words
        character(len=200):: words, temp_word
        integer :: i, j, ist

        ! remove space
        temp_word = trim(input_line)

        ! cout word
        num_words = 0
        ist = 1
        999 continue
        do i = ist, len_trim(temp_word)
            if (temp_word(i:i) /= ' ') then
                num_words = num_words + 1
                j = 0
                do while (i+j <= len_trim(temp_word) .and. temp_word(i+j:i+j) /= ' ')
                    j = j + 1
                end do
                words = temp_word(i:i+j-1)
                WRITE(*,*)trim(words)
                ist = i + j
                goto 999
            end if
        end do
  END SUBROUTINE split_string


