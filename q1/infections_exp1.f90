program infections_exp1 
    ! we will use the module in which we defined our model 
    use infections

    implicit none 

    ! we define single and double precision to easily switch between both 
    integer, parameter :: sp = kind(1.0) ! single precision 
    integer, parameter :: dp = kind(1.0d0)  ! double precision 
    ! following line to change for different precision 
    integer, parameter :: select_p= sp ! set select_p to the desired precision (sp or dp)
    real(select_p) :: r, d, result
    integer i

    ! we prompt the user for the 3 parameters of the model 
    print *,'Enter the number of initial infections'
    read *,i
    print *, 'Enter the number of days'
    read *,d
    print *, 'Enter the daily infection rate'
    read *, r

    ! catch problems 
    if (i < 0 .or. int(i) /= i) then 
        print *,"ERROR: the number of initial infections is invalid"
    elseif(d < 0 .or. int(d) /= d) then 
        print *, "ERROR: the number of days is invalid"
    elseif(abs(r) > 1) then 
        print *,"ERROR: rate is invalid"
    else 
        ! we call the function to compute the number of infections 
        result = predicted_infections(r,d,i)
        print *, "The number of predicted infections is ", result

    endif 
    
    
end program 