  PROGRAM post_monitor
  IMPLICIT NONE
  INTEGER::ni,nj,nk
  INTEGER::ist,ied,jst,jed,kst,ked
  REAL(8),ALLOCATABLE::x (:,:,:),y (:,:,:),z (:,:,:)
  REAL(8),ALLOCATABLE::xu(:,:,:),yu(:,:,:),zu(:,:,:)
  REAL(8),ALLOCATABLE::xv(:,:,:),yv(:,:,:),zv(:,:,:)
  REAL(8),ALLOCATABLE::xw(:,:,:),yw(:,:,:),zw(:,:,:)
  REAL(8),ALLOCATABLE::c(:,:,:),tpr(:,:,:),u(:,:,:),v(:,:,:),w(:,:,:)
  CHARACTER(len=32)::fname,fname_base
  INTEGER::nfile
  INTEGER::i,j,k
  CHARACTER(len=32)::ctmp

  fname_base="r4"

  ! read scalar grid & allocate
  fname = trim(fname_base) // "-scalar.mgrd"
  WRITE(*,*)"#Open file: ",trim(fname)
  OPEN(10,file=trim(fname),status='OLD')
    READ(10,*)ctmp,ist,ctmp,ied,ctmp,jst,ctmp,jed,ctmp,kst,ctmp,ked
    ni=ied-ist+1
    nj=jed-jst+1
    nk=ked-kst+1
    WRITE(*,*)"#Range: ni=",ni,"nj=",nj,"nk=",nk
    ALLOCATE(x(ni,nj,nk),y(ni,nj,nk),z(ni,nj,nk))
    ALLOCATE(xu(ni,nj,nk),yu(ni,nj,nk),zu(ni,nj,nk))
    ALLOCATE(xv(ni,nj,nk),yv(ni,nj,nk),zv(ni,nj,nk))
    ALLOCATE(xw(ni,nj,nk),yw(ni,nj,nk),zw(ni,nj,nk))
    ALLOCATE(c(ni,nj,nk),tpr(ni,nj,nk),u(ni,nj,nk),v(ni,nj,nk),w(ni,nj,nk))

    ! read grid
    READ(10,*)(((x(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((y(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((z(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
  CLOSE(10)

  ! read vector-U grid
  fname = trim(fname_base) // "-vector-U.mgrd"
  WRITE(*,*)"#Open file: ",trim(fname)
  OPEN(10,file=trim(fname),status='OLD')
    READ(10,*)
    ! read grid
    READ(10,*)(((xu(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((yu(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((zu(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
  CLOSE(10)

  ! read vector-V grid
  fname = trim(fname_base) // "-vector-V.mgrd"
  WRITE(*,*)"#Open file: ",trim(fname)
  OPEN(10,file=trim(fname),status='OLD')
    READ(10,*)
    ! read grid
    READ(10,*)(((xv(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((yv(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((zv(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
  CLOSE(10)

  ! read vector-W grid
  fname = trim(fname_base) // "-vector-W.mgrd"
  WRITE(*,*)"#Open file: ",trim(fname)
  OPEN(10,file=trim(fname),status='OLD')
    READ(10,*)
    ! read grid
    READ(10,*)(((xw(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((yw(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
    READ(10,*)(((zw(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
  CLOSE(10)

  fname = trim(fname_base) // "-tpr" // "_000011" // ".mval"
  OPEN(10,file=trim(fname),status='OLD')
    READ(10,*)
    ! x
    READ(10,*)(((tpr(i,j,k),k=kst,ked),j=jst,jed),i=ist,ied)
  CLOSE(10)

  OPEN(11,file='test-tpr.dat')
    WRITE(11,*)'VARIABLES="X" "Y" "Z" "tpr"'
    WRITE(11,*)'ZONE I=',ni,"J=",nj,"K=",nk
    WRITE(11,*)'DATAPACKING=POINT'
    DO k=kst,ked
    DO j=jst,jed
    DO i=ist,ied
      WRITE(11,'(4E12.4)')x(i,j,k),y(i,j,k),z(i,j,k),tpr(i,j,k)
    ENDDO
    ENDDO
    ENDDO
  CLOSE(11)

  STOP
  END
