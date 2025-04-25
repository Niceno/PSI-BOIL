subroutine cscs_read_procstatm()

 implicit none
 character(len=20) :: file
 integer :: size=0, resident=0, share=0, text=0, lib=0, data=0, dt=0
 integer:: ierr=0

 file='/proc/self/statm'
  open(unit=1,file=file,form='formatted',STATUS='OLD',ACTION='READ',IOSTAT=ierr)
  read(unit=1,FMT=*,IOSTAT=ierr) size, resident, share, text, lib, data, dt
 close(unit=1)

 if (ierr < 0) then
  write(*,*)'Problem reading /proc/self/statm'
 else
  write(*,'(a,i8)') '/proc/self/statm size = ', size
  write(*,'(a,i8)') '/proc/self/statm resi = ', resident
  write(*,'(a,i8)') '/proc/self/statm data = ', data
 end if
end subroutine cscs_read_procstatm

