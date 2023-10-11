implicit none
real :: a,b,c,fa,fb,fc,tol,g,t,x
integer :: N

tol = 1e-5
g = -32.17
t = 1
x = 1.7
a = -8
b = 0
fa = g * (cosh(a*t) - cos(a*t)) - (x/t) * 2 * a
fb = g * (cosh(b*t) - cos(b*t)) - (x/t) * 2 * b

do N = 1,50
    c = (a + b) / 2
    fc = g * (cosh(c*t) - cos(c*t)) - (x/t) * 2 * c
    print*,'a = ',a,'b = ',b,'c = ',c
    if (fc == 0 .OR. (fc > 0 .AND. fc <= tol) .OR. (fc < 0 .AND. fc >= (-1)*tol)) then
        print*,'approximately zero at w = ',c, 'f(w) = ',fc
        stop
    end if
    if (fc > 0 .AND. fa > 0) then
        a = c
        fa = g * (cosh(a*t) - cos(a*t)) - (x/t) * 2 * a
    else if (fc < 0 .AND. fa < 0) then
        a = c
        fa = g * (cosh(a*t) - cos(a*t)) - (x/t) * 2 * a
    else
        b = c
        fb = g * (cosh(b*t) - cos(b*t)) - (x/t) * 2 * b
    end if
end do

stop
end




