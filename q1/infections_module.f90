module infections 
    implicit none 
    contains 
    ! function: predicted infections 
    ! computes the number of infections according to 3 parameters: 
    ! r: the infection rate 
    ! d: the number of days 
    ! i0: initial number of infections 
    function predicted_infections(r,d,i0) result(id)
        ! we define single and double precision to easily switch between both 
        integer, parameter :: sp = kind(1.0) ! single precision 
        integer, parameter :: dp = kind(1.0d0)  ! double precision 
        ! following line to change for different precision 
        integer, parameter :: select_p=sp ! set select_p to the desired precision (sp or dp)
        real(select_p) :: r, d, id
        integer i0 
        id = (1+r)**d * i0 ! computing infections according to the model 
    end function 

end module 