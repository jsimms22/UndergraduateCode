program num4_2
implicit none
integer :: a,b,n

a = 0
b = 2
n = 10000
call trapezoidal(a,b,n)
call simpson(a,b,n)
call romberg(a,b,n)

stop
end

subroutine trapezoidal(a,b,n)
implicit none
integer :: i,a,b,n
real*8 :: h,t1,t2,integral,j

h = ( dble(b) - dble(a) ) / dble(n)
t1 = (1.0 - dble(a)**3) + ((1.0 - dble(b))**2)
t2 = 0.0

!print*,'first loop'
do i = 1,(n-1)
    j = dble(a) + dble(i)*h
    if (j <= 1) then
        t2 = t2 + (1.0 - j**3)
    else
        t2 = t2 + ((1.0 - j)**2)
    endif
    !print*,i,j,t2
enddo

integral =  (h/2.0)*(t1 + 2.0*t2)

print*,'trapezoidal approx = ',integral

return
end


subroutine simpson(a,b,n)
implicit none
integer :: i,a,b,n
real*8 :: h,t1,t2,t3,integral,j

h = ( real(b) - real(a) ) / real(n)
t1 = (1.0 - real(a)**3) + ((1.0 - real(b))**2)
t2 = 0.0
t3 = 0.0

!print*,'first loop'
do i = 1,(n-1),2
    j = real(a) + real(i)*h
    if (j <= 1) then
        t2 = t2 + (1.0 - j**3)
    else
        t2 = t2 + ((1.0 - j)**2)
    endif
    !print*,i,j,t2
enddo
!print*,'2nd loop'
do i = 2,n,2
    j = real(a) + real(i)*h
    if (j <= 1) then
        t3 = t3 + (1.0 - j**3)
    else
        t3 = t3 + ((1.0 - j)**2)
    endif
    !print*,i,j,t3
enddo

integral =  (h/3.0)*(t1 + 2.0*t2 + 4.0*t3)

print*,'simpson approx = ',integral

return
end

subroutine romberg(a,b,n)
implicit none
integer :: i,j,k,a,b,n
real*8 :: h,kk
real*8, dimension(2,n) :: R

h = real(b) - real(a)
R(1,1) = (h/2.0)*((1.0 - real(a)**3) + ((1.0 - real(b))**2))
!print*,'R(1,1) = ',R(1,1)
do i = 2,n
    R(2,1) = 0.5*(R(1,1))
    do k = 1,2**(i-2)
        kk = real(a) + (real(k) - 0.5)*h
        if (kk <= 1) then
            R(2,1) = R(2,1) + (h/2.0)*(1.0 - kk**3)
        else
            R(2,1) = R(2,1) + (h/2.0)*((1.0 - kk)**2)
        endif
    enddo
    !print*,'R(2,1) = ',R(2,1)
    do j = 2,i
        R(2,j) = R(2,j-1) + (R(2,j-1) - R(1,j-1))/(4.0**(j-1) - 1)
        !print*,'R(2,j) = ',R(2,j)
    enddo
    h = h/2.0
    do j = 1,i
        if (sqrt((R(2,i) - R(1,i-1))**2) <= 1e-8 ) then
            print*,'exiting romberg subroutine early, answer within tolerance'
            print*,'romberg approx = ',R(2,i)
            return
        endif
        R(1,j) = R(2,j)
    enddo
enddo

print*,'romberg approx = ',R(2,n)
    
return
end
