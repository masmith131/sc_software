program test_solver
    use solver
    implicit none 
    
    ! make a matrix A 
    real(wp), dimension(3,3) :: A
    real(wp), dimension(3) ::b

    real(wp), dimension(5,5) :: M
    real(wp), dimension(5) :: eigenvalues
    integer i, j

    ! A = reshape([1.0,1.0, 2.0, 2.0, 1.0, 3.0, 1.0, 4.0,2.0], [3, 3])
    A = reshape([1.0,2.0,1.0,1.0,1.0,4.0,2.0,3.0,2.0], [3,3])
    b = [8.0,13.0,11.0]
    call solve(A,b)
    print *,b

    ! Define a 5x5 matrix A
    M = reshape([ &
        1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, &
        2.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, &
        3.0_wp, 3.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, &
        4.0_wp, 4.0_wp, 4.0_wp, 4.0_wp, 5.0_wp, &
        5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp], &
        [5, 5])

    ! Print the matrix
    do i = 1,5
        do j = 1,5
            write(*, '("(",F6.2,",",A)', advance='no') M(i,j)
        end do
        write(*,*)  ! New line for each row
    end do

    ! Call the get_eigenvalues subroutine
    call get_eigenvalues(M, eigenvalues)

    ! Output the eigenvalues
    print *, 'The eigenvalues of matrix M are: ', eigenvalues



end program 