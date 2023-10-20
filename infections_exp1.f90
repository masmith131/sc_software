program infections_exp1 
    use infections

    implicit none 
    integer, parameter :: sp = kind(1.0) ! single precision 
    integer, parameter :: dp = kind(1.0d0)  ! double precision 
    integer, parameter :: select_p= dp ! set select_p to the desired precision 
    real(select_p) :: r, d, result
    integer i

    print *,'Enter the number of initial infections'
    read *,i
    print *, 'Enter the number of days'
    read *,d
    print *, 'Enter the daily infection rate'
    read *, r
    
    result = predicted_infections(r,d,i)
    print *, "The number of predicted infections is ", result
    !if (i > 1e9) then 
    !    print *, 'The number you entered is too big :(' 
    !elseif (i < 0) then 
    !    print *, 'You entered a negative number :('
    !endif
    
end program 