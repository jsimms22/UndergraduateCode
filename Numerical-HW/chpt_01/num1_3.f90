implicit none
real :: p,temp,val
integer :: ii

temp = 0.0
p = 0.1
val = 0.0

do
    temp = p
    temp = 1.0 + p
    p = p / 1.01
    val = 1.0 + p
    !!print*,val
    if (val == 1.0) then
        print*,'Closest estimated real value to 1: ',temp
        print*,'Processor: Intel Core i7-5500U CPU @ 2.40GHz'
        print*,'I am using a virtual box, so only half my cores are'
        print*,'used'
        print*,'Compiler: gfortran'
        stop
    end if
end do

stop
end
