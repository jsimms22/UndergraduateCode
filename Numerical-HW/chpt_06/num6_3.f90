program num6_3
implicit none
integer :: steps
real*8 :: initpos,initvel,PI,a,b

PI = 4.0*ATAN(1.0d0)
steps = 100
a = 0
b = 24.0d0*PI
initpos = 1.0d0
initvel = 0.0d0

call modeuler(a,b,steps,initpos,initvel)

!call ralston(a,b,steps,init)

!call rk4(a,b,steps,init)

stop
end

real*8 function accel(u,theta)
implicit none
real*8 :: theta,u
accel = 0.0d0
accel = 1.0d0 - u  
return
end

real*8 function x(u,theta)
implicit none
real*8 :: theta,u
x = (1.0d0/u)*COS(theta)
return
end

real*8 function y(u,theta)
implicit none
real*8 :: theta,u
y = (1.0d0/u)*SIN(theta)
return
end

subroutine modeuler(a,b,steps,initpos,initvel)
implicit none
integer :: steps,ii
real*8 :: initpos,initvel,tstep,theta,accel,a,b,x,y
real*8, dimension(steps+1) :: u,v

open(unit=100,file="orbit.dat",form="formatted",status="replace",action="write")

tstep = (dble(b) - dble(a)) / dble(steps)
theta = dble(a)
u(1) = initpos
v(1) = initvel

write(unit=100,fmt=*) x(u(1),theta),y(u(1),theta) !theta,u(1),v(1),accel(u(1),theta)

do ii = 2,steps+1
    v(ii) = v(ii-1) + tstep * accel(u(ii-1),theta)
    u(ii) = u(ii-1) + tstep * v(ii)
    v(ii) = v(ii-1) + tstep/2.0d0 * (accel(u(ii-1),theta) + accel(u(ii),theta+tstep))
    u(ii) = u(ii-1) + tstep/2.0d0 * (v(ii-1) + v(ii))
    theta = a + (ii-1)*tstep
    write(unit=100,fmt=*) x(u(ii),theta),y(u(ii),theta) !theta,u(ii),v(ii),accel(u(ii),theta)
enddo

close(unit=100)

return
end

