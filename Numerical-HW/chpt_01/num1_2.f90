implicit none
real :: small, temp

temp = 0.0
small = 0.1

do
    temp = small
    small = small / 10.0
    !!print*,small
    if (small == 0.0) then
        print*,'Smallest estimated real value: ',temp
        stop
    end if
end do

stop
end
