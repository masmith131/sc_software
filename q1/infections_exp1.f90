program infections_exp1 
    ! we will use the module in which we defined our model 
    use infections

    implicit none 

    ! use precision select_p defines in the module infections
    real(select_p) :: r, d, result ! rate, number of days, final result 
    integer i ! number of initial infections

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