program infections_exp2
    use infections

    implicit none 
    ! we define single and double precision to easily switch between both 
    integer, parameter :: sp = kind(1.0) ! single precision 
    integer, parameter :: dp = kind(1.0d0)  ! double precision 
    ! following line to change for different precision 
    integer, parameter :: select_p= dp ! set select_p to the desired precision (sp or dp)
    real(select_p) :: r, d, result
    integer i 

    ! define here the rate r, the number of days d and the initial number of infections i 
    r = 0.01 ! change this value to r = 0.009765625 for Q1.3
    d = 1500.5
    i = 100
    
    ! we call the function to compute the number of infections 
    result = predicted_infections(r,d,i)
    print *, "The result of this experiment is ", result
    
end program 