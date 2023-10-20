module infections 
    implicit none 
    contains 
    function predicted_infections(r,d,i0) result(id)
        integer, parameter :: sp = kind(1.0) ! single precision 
        integer, parameter :: dp = kind(1.0d0)  ! double precision 
        integer, parameter :: select_p= dp ! set select_p to the desired precision 
        real(select_p) :: r, d, id
        integer i0 
        id = (1+r)**d * i0
    end function 

end module 