program test_solver
    use solver
    implicit none 
    
    ! make a matrix A 
    real, dimension(3,3) :: A
    real, dimension(3) ::b

    ! A = reshape([1.0,1.0, 2.0, 2.0, 1.0, 3.0, 1.0, 4.0,2.0], [3, 3])
    A = reshape([1.0,2.0,1.0,1.0,1.0,4.0,2.0,3.0,2.0], [3,3])
    b = [8.0,13.0,11.0]
    call solve(A,b)
    print *,b


end program 