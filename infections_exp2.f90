program infections_exp2
    use infections

    implicit none 
    integer, parameter :: sp = kind(1.0) ! single precision 
    integer, parameter :: dp = kind(1.0d0)  ! double precision 
    integer, parameter :: select_p= dp ! set select_p to the desired precision 
    real(select_p) :: r, d, result
    integer i 

    r = 0.01
    d = 1500.5
    i = 100
    
    result = predicted_infections(r,d,i)
    print *, "The result of this experiment is ", result
    
end program 