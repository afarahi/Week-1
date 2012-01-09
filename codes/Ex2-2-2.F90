program Ex2_2_2

real    :: x0 , x1 , x2 , x
integer :: i 

x0 = 1.0d0
x1 = 4.0d0

do i = 2 , 15
   x2 = 5.0d0 * x1 - 4.0d0 * x0
   x0 = x1
   x1 = x2
   write (*,*) x2
end do

x = (4.0d0)**15

write (*,*) 'Relative error is : ' , ( abs(x2 - x) / x ) 

write (*,*) 'Absolute error is : ' , abs(x2 - x)

end program Ex2_2_2

