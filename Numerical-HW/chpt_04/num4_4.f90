program num4_3
implicit none
integer :: m,n
real*8 :: a,b,c,d,bigM,bigMy,bigMx,xbar,ybar

m = 14
n = 12
a = 0
b = 3.14159 / 4.0
print*,'m = ',m,'n = ',n
print*,'M = ',bigM(a,b,m,n)
print*,'Mx = ',bigMx(a,b,m,n)
print*,'My = ',bigMy(a,b,m,n)
xbar = bigMx(a,b,m,n)/bigM(a,b,m,n)
ybar = bigMy(a,b,m,n)/bigM(a,b,m,n)

print*,'xbar = ',xbar,'ybar = ',ybar
print*,''
m = 100
n = 100
a = 0
b = 3.14159 / 4.0
print*,'m = ',m,'n = ',n
print*,'M = ',bigM(a,b,m,n)
print*,'Mx = ',bigMx(a,b,m,n)
print*,'My = ',bigMy(a,b,m,n)
xbar = bigMx(a,b,m,n)/bigM(a,b,m,n)
ybar = bigMy(a,b,m,n)/bigM(a,b,m,n)

print*,'xbar = ',xbar,'ybar = ',ybar


stop
end

real*8 function bigM(a,b,m,n)
implicit none
integer :: m,n,ii,jj
real*8 :: a,b,c,d,h,hx
real*8 :: J,J1,J2,J3,y
real*8 :: x,K1,K2,K3,Q,L

h = (b - a) / dble(n)
J1 = 0.0
J2 = 0.0
J3 = 0.0

do ii = 0,n
    x = a + dble(ii)*h
    hx = (cos(x) - sin(x)) / dble(m)
    K1 = sqrt(x)*exp(x**2 + sin(x)**2) + sqrt(x)*exp(x**2 + cos(x)**2)
    K2 = 0.0
    K3 = 0.0
    do jj = 1,m-1
        y = sin(x) + jj*hx
        Q = sqrt(x) * exp(x**2 + y**2)
        if ( MOD(jj,2) == 0 ) then
            K2 = K2 + Q
        else
            K3 = K3 + Q
        endif
        L = (K1 + 2.0*K2 + 4.0*K3) * (hx/3.0)
        if (ii == 0 .OR. ii == n) then
            J1 = J1 + L
        elseif ( MOD(ii,2) == 0 ) then
            J2 = J2 + L
        else
            J3 = J3 + L
        endif
    enddo
enddo

bigM = h * (J1 + 2.0*J2 + 4.0*J3) / 3.0

return
end

real*8 function bigMy(a,b,m,n)
implicit none
integer :: m,n,ii,jj
real*8 :: a,b,c,d,h,hx
real*8 :: J,J1,J2,J3,y
real*8 :: x,K1,K2,K3,Q,L

h = (b - a) / dble(n)
J1 = 0.0
J2 = 0.0
J3 = 0.0

do ii = 0,n
    x = a + dble(ii)*h
    hx = (cos(x) - sin(x)) / dble(m)
    K1 = sin(x) * sqrt(x)*exp(x**2 + sin(x)**2) + cos(x) * sqrt(x)*exp(x**2 + cos(x)**2)
    K2 = 0.0
    K3 = 0.0
    do jj = 1,m-1
        y = sin(x) + jj*hx
        Q = y * sqrt(x) * exp(x**2 + y**2)
        if ( MOD(jj,2) == 0 ) then
            K2 = K2 + Q
        else
            K3 = K3 + Q
        endif
        L = (K1 + 2.0*K2 + 4.0*K3) * (hx/3.0)
        if (ii == 0 .OR. ii == n) then
            J1 = J1 + L
        elseif ( MOD(ii,2) == 0 ) then
            J2 = J2 + L
        else
            J3 = J3 + L
        endif
    enddo
enddo

bigMy = h * (J1 + 2.0*J2 + 4.0*J3) / 3.0

return
end

real*8 function bigMx(a,b,m,n)
implicit none
integer :: m,n,ii,jj
real*8 :: a,b,c,d,h,hx
real*8 :: J,J1,J2,J3,y
real*8 :: x,K1,K2,K3,Q,L

h = (b - a) / dble(n)
J1 = 0.0
J2 = 0.0
J3 = 0.0

do ii = 0,n
    x = a + dble(ii)*h
    hx = (cos(x) - sin(x)) / dble(m)
    K1 = x * sqrt(x)*exp(x**2 + sin(x)**2) + x * sqrt(x)*exp(x**2 + cos(x)**2)
    K2 = 0.0
    K3 = 0.0
    do jj = 1,m-1
        y = sin(x) + jj*hx
        Q = x * sqrt(x) * exp(x**2 + y**2)
        if ( MOD(jj,2) == 0 ) then
            K2 = K2 + Q
        else
            K3 = K3 + Q
        endif
        L = (K1 + 2.0*K2 + 4.0*K3) * (hx/3.0)
        if (ii == 0 .OR. ii == n) then
            J1 = J1 + L
        elseif ( MOD(ii,2) == 0 ) then
            J2 = J2 + L
        else
            J3 = J3 + L
        endif
    enddo
enddo

bigMx = h * (J1 + 2.0*J2 + 4.0*J3) / 3.0

return
end
