program num5_1
implicit none
integer :: steps,ii
real*8 :: x,y,PI,sine

open(unit=7,file="sine.dat",form="formatted",status="replace",action="write")

PI = 3.14
steps = 200

do ii = 0,steps
    x = dble(ii) * (1.0d0/100.0d0) * PI
    y = sine(x)
    write(unit=7,fmt=*) y
enddo

close(unit=7)

stop
end

real*8 function sine(x)
implicit none
real*8 :: x

sine = sin(x)
print*,sine

return
end
