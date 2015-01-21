program dipfft
implicit none
!
complex*16 ::s,Ew 
complex*16 :: dept_dipx(1:600201)
real*8::Ex(1:600201)
real*8::nmomega,omega,omega_start,omega_end,domega,dep,pi 
real*8::hq,anf,hq1,temp1,temp2
real*8:: dtt,tt_min,tt_max
real*8:: tt(1:600201)
real*8:: time,t_dipx
complex*16:: f_dipx
integer*4:: i,j,k,line_output
integer*4:: line_input
character*30:: file_output,file_input,E_input
namelist /filename/ file_input,E_input,file_output
namelist /freq/omega_start,omega_end,line_input,line_output,dep

!
rewind 5
read(5,filename,end=101)
101 continue
PRINT *,'input file name',file_input
PRINT *,'input E file name',E_input
PRINT *,'output file name',file_output


rewind 5
read(5,freq,end=102)
102 continue 
PRINT *,'starting freq. in eV',omega_start
PRINT *,'end freq. in eV',omega_end
PRINT *,'number of lines in the file',line_input
PRINT *,'number of output points',line_output
PRINT *,'dephasing',dep
!
hq=0.65822d0
hq1=hq
anf=0.00000d0
pi=3.141592653d0
!
dep=-1d0*dep/hq
!



!
OPEN(70,FILE=file_input,STATUS='OLD')
do i=1,line_input
 read(70,*)time,t_dipx
  tt(i)=time
  dept_dipx(i)=(1.,0.)*t_dipx*dexp(dep*time)
enddo
close(70)

open(71,FILE=E_input,STATUS='OLD')
do i=1,line_input
read(71,*)temp1,temp2
Ex(i)=(1.,0.)*temp2
enddo
close(71)

	
!
OPEN(72,FILE=file_output,STATUS='Replace')
tt_min=tt(1)
tt_max=tt(line_input)
dtt=tt(2)-tt(1)
domega=(omega_end-omega_start)/(line_output-1.d0)
omega=omega_start-domega
!
do i=1,line_output
  omega=omega+domega
  s=(0.d0,0.d0)
  Ew=(0.d0,0.d0)
  do j=1,line_input-1
    dtt=tt(j+1)-tt(j)
    s=s+exp(dcmplx(0.d0,omega*tt(j)/hq))*dept_dipx(j)*dtt
    Ew=Ew+exp(dcmplx(0.d0,omega*tt(j)/hq))*Ex(j)*dtt
  enddo
  omega=omega
         f_dipx=s/Ew/(27.2113986955332d0,0.d0)                         
 
 WRITE(72,'(5E24.14)') omega,real(f_dipx),aimag(f_dipx),abs(f_dipx),omega*aimag(f_dipx)*2.0/3.0/pi
enddo
close(72)
end program dipfft
