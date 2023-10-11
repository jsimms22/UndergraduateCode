program testdivdiff
implicit none
integer :: i,n
real*8, allocatable, dimension(:) :: xdata,ydata,a
real*8 :: x,p

n = 5
allocate(xdata(n),ydata(n),a(n))

xdata = (/ 1.0d0,2.0d0,3.0d0,4.0d0,5.0d0 /)
ydata = (/ -4.0d0,-11.0d0,-14.0d0,-7.0d0,16.0d0 /)

call divdiff(n,xdata,ydata,a)

do i = 1,n
    x = xdata(i)
    print*,xdata(i),ydata(i),ydata(i) - p(n,xdata,a,x)
enddo

deallocate(xdata,ydata,a)

stop
end

!Subroutine divdiff uses Newton's divided difference
!algorthim to calculate the coefficients of an
!interpolation polynominal in the Newton form.
!
!   n is number of data points
!   xdata is the values of the data's x components
!   ydata is the values of the data's y components
!   a is a 1-dimensional array to store the coefficients
!       before returning to the main program
!
!Divdiff will print the 1-d array before returning and
!ending the subroutine

subroutine divdiff(n,xdata,ydata,a)
implicit none
integer n,ii,jj
real*8, dimension(n) :: a,xdata,ydata
real*8, allocatable, dimension(:,:) :: z

allocate(z(n,n))
z = 0
a = 0
do ii = 1,n
    z(ii,1) = ydata(ii)
enddo

a(1) = z(1,1)

do jj = 2,n
    do ii = 1,n
        z(ii,jj) = (z(ii+1,jj-1) - z(ii,jj-1)) / (xdata(ii+jj-1) - xdata(ii))
    enddo
    a(jj) = z(1,jj)
enddo

print*,'Coefficients from divdiff'
do ii = 1,n
    print*,a(ii)
enddo
    print*,

deallocate(z)

return
end

real*8 function p(n,xdata,a,x)
implicit none
integer :: ii,n
real*8 :: x,xdata(n),a(n)

p = a(n)

do ii = n,2,-1
    p = p*(x - xdata(ii-1)) + a(ii-1)
enddo

return
end
