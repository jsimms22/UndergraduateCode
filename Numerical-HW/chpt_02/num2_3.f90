implicit none
real :: a,b,c,fa,fb,fc,tol,P
integer :: N
P = .5
tol = 1e-3
a = 0
b = 1
fa = ((1+a)/2)*(a/(1-a+a**2))**21 - P
fb = ((1+b)/2)*(b/(1-b+b**2))**21 - P
 
do N = 1,50
    c = (a + b) / 2
    fc = ((1+c)/2)*(c/(1-c+c**2))**21 - P
    print*,'a = ',a,'b = ',b,'c = ',c
    if (fc == 0 .OR. (fc > 0 .AND. fc <= tol) .OR. (fc < 0 .AND. fc >= (-1)*tol)) then
        print*,'minimum p value = ',c
        stop
    end if
    if (fc > 0 .AND. fa > 0) then
        a = c
        fa = ((1+a)/2)*(a/(1-a+a**2))**21 - P
    else if (fc < 0 .AND. fa < 0) then
        a = c
        fa = ((1+a)/2)*(a/(1-a+a**2))**21 - P
    else
        b = c
        fb = ((1+b)/2)*(b/(1-b+b**2))**21 - P
    end if
end do

stop
end

