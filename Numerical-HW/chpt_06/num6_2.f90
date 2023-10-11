program num6_2
implicit none
integer :: steps
real*8 :: init,a,b

steps = 100
a = 0.0d0
b = 0.20
init = 1.0d0

call modeuler(a,b,steps,init)

call ralston(a,b,steps,init)

stop
end

real*8 function yprime(y,t)
implicit none
real*8 :: t,y
yprime = 0.0d0
yprime = exp(y)
return
end

subroutine modeuler(a,b,steps,init)
implicit none
integer :: steps,ii
real*8 :: init,tstep,t,yprime,a,b
real*8, dimension(steps+1) :: y

open(unit=7,file="modeuler2.dat",form="formatted",status="replace",action="write")

tstep = (dble(b) - dble(a)) / dble(steps)
t = dble(a)
y(1) = dble(init)

write(unit=7,fmt=*) t,y(1)

do ii = 2,steps+1
    y(ii) = y(ii-1) + tstep * yprime(y(ii-1),t) 
    y(ii) = y(ii-1) + tstep/2.0d0 * (yprime(y(ii-1),t) + yprime(y(ii),t+tstep))
    t = a + (ii-1)*tstep
    write(unit=7,fmt=*) t,y(ii)
enddo

close(unit=7)

return
end


subroutine ralston(a,b,steps,init)
implicit none
integer :: steps,ii
real*8 :: init,tstep,k1,k2,t,yprime,a,b
real*8, dimension(steps+1) :: y

open(unit=8,file="ralston2.dat",form="formatted",status="replace",action="write")

tstep = (dble(b) - dble(a)) / dble(steps)
t = dble(a)
y(1) = dble(init)

write(unit=8,fmt=*) t,y(1)

do ii = 2,steps+1
    k1 = yprime(y(ii-1),t)
    k2 = yprime((y(ii-1) + (2.0d0/3.0d0) * tstep * k1),(t + tstep*(2.0d0/3.0d0)))
    y(ii) = y(ii-1) + tstep * ((1.0d0/4.0d0) * k1 + (3.0d0/4.0d0) * k2)
    t = a + (ii-1) * tstep
    write(unit=8,fmt=*) t,y(ii)
enddo

close(unit=8)

return
end


