program num5_3
implicit none
integer :: steps,a,b
real*8 :: init

steps = 100
a = 0
b = 2
init = 200.0d0

call euler(a,b,steps,init)

call ralston(a,b,steps,init)

call rk4(a,b,steps,init)

stop
end

subroutine euler(a,b,steps,init)
implicit none
integer :: a,b,steps,ii
real*8 :: init,tstep,t
real*8, dimension(steps+1) :: y

open(unit=7,file="euler.dat",form="formatted",status="replace",action="write")

tstep = (dble(b) - dble(a)) / dble(steps)
t = dble(a)
y(1) = dble(init)

write(unit=7,fmt=*) t,y(1)

do ii = 2,steps+1
    y(ii) = y(ii-1) + (tstep * (t * exp(3.0d0 * t) - 2.0d0 * y(ii-1)))
    t = a + (ii-1)*tstep
    write(unit=7,fmt=*) t,y(ii)
enddo

close(unit=7)

return
end

subroutine ralston(a,b,steps,init)
implicit none
integer :: a,b,steps,ii
real*8 :: init,tstep,k1,k2,t
real*8, dimension(steps+1) :: y

open(unit=8,file="ralston.dat",form="formatted",status="replace",action="write")

tstep = (dble(b) - dble(a)) / dble(steps)
t = dble(a)
y(1) = dble(init)

write(unit=8,fmt=*) t,y(1)

do ii = 2,steps+1
    k1 = t * exp(3.0d0 * t) - 2.0d0 * y(ii-1)
    k2 = (t + tstep*(2.0d0/3.0d0)) * exp(3.0d0 * (t + tstep*(2.0d0/3.0d0))) - 2.0d0 * (y(ii-1) + (2.0d0/3.0d0) * tstep * k1)
    y(ii) = y(ii-1) + tstep * ((1.0d0/4.0d0) * k1 + (3.0d0/4.0d0) * k2)
    t = a + (ii-1) * tstep
    write(unit=8,fmt=*) t,y(ii)
enddo

close(unit=8)

return
end

subroutine rk4(a,b,steps,init)
implicit none
integer :: a,b,steps,ii
real*8 :: init,tstep,k1,k2,k3,k4,t
real*8, dimension(steps+1) :: y

open(unit=9,file="rk4.dat",form="formatted",status="replace",action="write")

tstep = (dble(b) - dble(a)) / dble(steps)
t = dble(a)
y(1) = dble(init)

write(unit=9,fmt=*) t,y(1)

do ii = 2,steps+1
    k1 = tstep * (t * exp(3.0d0 * t) - 2.0d0 * y(ii-1))
    k2 = tstep * ((t + tstep/2.0d0) * exp(3.0d0 * (t + tstep/2.0d0)) - 2.0d0 * (y(ii-1) + k1/2.0d0))
    k3 = tstep * ((t + tstep/2.0d0) * exp(3.0d0 * (t + tstep/2.0d0)) - 2.0d0 * (y(ii-1) + k2/2.0d0))
    k4 = tstep * ((t + tstep) * exp(3.0d0 * (t + tstep)) - 2.0d0 * (y(ii-1) + k3))
    y(ii) = y(ii-1) + (1.0d0/6.0d0)*(k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)
    t = a + (ii-1) * tstep
    write(unit=9,fmt=*) t,y(ii)
enddo

close(unit=9)

return
end
