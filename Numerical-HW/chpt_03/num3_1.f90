implicit none
integer :: i, n = 5
real*8, dimension(5) :: x,y,a

do i = 1,n
    x(i) = i
    y(i) = i*i*i
enddo

call divdiff(n,x,y,a)

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
z = 0;
a = 0;
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

do ii = 1,n
    print*,a(ii)
enddo

deallocate(z)

return
end
