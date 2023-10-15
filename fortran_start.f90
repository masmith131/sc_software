program question1 
    implicit none 
    integer i 
    print *,'Enter a number'
    read *,i
    if (i > 1e9) then 
        print *, 'The number you entered is too big :(' 
    elseif (i < 0) then 
        print *, 'You entered a negative number :('
    endif
end program