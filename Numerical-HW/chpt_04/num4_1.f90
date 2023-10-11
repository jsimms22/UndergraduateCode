program fourone
implicit none
integer :: no
real*8 :: po,tol,a
no = 100
po = .05
a = 0.0
tol = 1e-5

call newtonNsimpson(po,tol,no,a)
call newtonNtrapezoidal(po,tol,no,a)

stop
end

real*8 function simpson(a,b,n)
implicit none
integer :: i,n
real*8 :: h,t1,t2,t3,integral,j,a,b

h = ( b - a ) / dble(n)
t1 = (1.0/sqrt(2.0 * 3.14159))*exp((-1.0*a**2)/2.0) + (1.0/sqrt(2.0 * 3.14159))*exp((-1.0*b**2)/2.0)
t2 = 0.0
t3 = 0.0

!print*,'first loop'
do i = 1,(n-1),2
    j = a + dble(i)*h
    t2 = (1.0/sqrt(2.0 * 3.14159))*exp(j**2/2.0)
    !print*,i,j,t2
enddo
!print*,'2nd loop'
do i = 2,n,2
    j = a + dble(i)*h
    t3 = (1.0/sqrt(2.0 * 3.14159))*exp(j**2/2.0)
    !print*,i,j,t3
enddo

simpson = (h/3.0)*(t1 + 2.0*t2 + 4.0*t3)
print*,simpson

return
end

real*8 function trapezoidal(a,b,n)
implicit none
integer :: i,n
real*8 :: h,t1,t2,j,a,b

h = ( b - a ) / dble(n)
t1 = (1.0/sqrt(2.0 * 3.14159))*exp((-1.0*a**2)/2.0) + (1.0/sqrt(2.0 * 3.14159))*exp((-1.0*b**2)/2.0)
t2 = 0.0

!print*,'first loop'
do i = 1,(n-1)
    j = a + dble(i)*h
    t2 = (1.0/sqrt(2.0 * 3.14159))*exp(j**2/2.0)
    !print*,i,j,t2
enddo

trapezoidal = (h/2.0)*(t1 + 2.0*t2)
print*,trapezoidal

return
end

subroutine newtonNsimpson(po,tol,no,a)
implicit none
integer :: i, no
real*8 :: po,p,tol,a,simpson

i = 1
do while (i <= no)
    p = po - (simpson(a,po,12)-0.45)/((1.0/sqrt(2.0 * 3.14159))*exp(po**2/2.0))
    if (sqrt((p-po)**2) <= tol) then
        print*,'exiting early, using simpsons rule p = ',p
        return
    endif
    i = i + 1
    po = p
enddo

print*,'method failed after ',no,' iterations'

return
end

subroutine newtonNtrapezoidal(po,tol,no,a)
implicit none
integer :: i, no
real*8 :: po,p,tol,a,trapezoidal

i = 1
do while (i <= no)
    p = po - (trapezoidal(a,po,4)-0.45)/((1.0/sqrt(2.0 * 3.14159))*exp(po**2/2.0))
    if (sqrt((p-po)**2) <= tol) then
        print*,'exiting early, using trapezoidal p = ',p
        return
    endif
    i = i + 1
    po = p
enddo

print*,'method failed after ',no,' iterations'

return
end
