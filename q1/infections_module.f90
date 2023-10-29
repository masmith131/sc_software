module infections 
    implicit none  

    ! we define single and double precision to easily switch between both 
    integer, parameter, public :: sp = kind(1.0) ! single precision 
    integer, parameter, public :: dp = kind(1.0d0)  ! double precision 
    ! following line to change for different precision 
    integer, parameter, public :: select_p=sp ! set select_p to the desired precision (sp or dp)

    contains
    ! =======================================================================================
    ! function: predicted_infections 
    ! computes the number of infections according to 3 parameters: 
    ! r: the infection rate 
    ! d: the number of days 
    ! i0: initial number of infections 
    ! ========================================================================================
    function predicted_infections(r,d,i0) result(id)
        real(select_p), intent(in) :: r, d
        integer, intent(in) :: i0 
        real(select_p) :: id 
        id = (1+r)**d*i0 ! computing infections according to the model 
    end function 
end module