program num4_3
implicit none
integer :: m,n
real*8 :: a,b,c,d

m = 2
n = 2
a = 0
b = 1
c = 0
d = 1

call simpson(a,b,c,d,m,n)

stop
end

subroutine simpson(a,b,c,d,m,n)
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
    hx = (d - c) / dble(m)
    K1 = x*c + x*d
    K2 = 0.0
    K3 = 0.0
    do jj = 1,m-1
        y = c + jj*hx
        Q = x*y
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

J = h * (J1 + 2.0*J2 + 4.0*J3) / 3.0

print*,'J = ',J

return
end
