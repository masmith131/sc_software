program infections_exp2
    use infections

    implicit none 
    ! use precision select_p defines in the module infections
    real(select_p) :: r_exp2, r_exp3, d, result
    integer i 

    ! define here the rate r, the number of days d and the initial number of infections i 
    r_exp2 = 0.01 ! rate used in experience 2 
    r_exp3 = 0.009765625 ! rate used in experience 3 
    d = 1500.5 ! number of days 
    i = 100 ! initial number of infections
    
    ! we call the function to compute the number of infections 
    result = predicted_infections(r_exp2,d,i) ! use r_exp2 or r_exp3
    ! print '(es20.10, 1x, a)', result, "is the result of number of infections"
    print *, result 
    
end program 