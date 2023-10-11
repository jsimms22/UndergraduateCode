implicit none
integer :: maxInt, temp

temp = 1;
maxInt = 1;

do
    temp = maxInt
    maxInt = maxInt + 1
    if (temp > maxInt) then
        print*,'Largest integer value is: ',temp
        stop
    end if
end do

stop
end
