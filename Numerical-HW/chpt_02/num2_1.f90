implicit none

real :: a,b,c,fa,fb,fc,tol
integer :: N
!f(x) = x**3 - 30
tol = 1e-5
a = 0
b = 7
fa = a**3 - 30
fb = b**3 - 30

do N = 1,50
    c = (a + b) / 2
    fc = (c)**3 - 30
    print*,'a = ',a,'b = ',b,'c = ',c
    !print*,'f(a) = ',fa,'f(b) = ',fb,'f(c) = ',fc
    if (fc == 0 .OR. (fc > 0 .AND. fc <= tol) .OR. (fc < 0 .AND. fc >= (-1)*tol)) then
        print*,'approximately zero at x = ',c, 'f(c) = ',fc
        stop
    end if
    if (fc > 0 .AND. fa > 0) then
        a = c
        fa = a**3 - 30
    else if (fc < 0 .AND. fa < 0) then
        a = c
        fa = a**3 - 30
    else
        b = c
        fb = b**3 - 30
    end if
end do

stop
end

