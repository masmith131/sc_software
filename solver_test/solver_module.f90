module solver
    implicit none

    integer, parameter :: s_p = kind(1.0) ! single precision 
    integer, parameter :: d_p = kind(1.0d0)  ! double precision 
    ! following line to change for different precision 
    integer, parameter, public :: wp=s_p ! set select_p to the desired precision (sp or dp)
    
    contains

    !---------------------------------------------------------------------------------------------
    ! solve(A,bx) solves the linear system Ax = b
    ! bx contains b at function call 
    ! bx contains x after function was exectuted
    ! it uses the Lapack subroutine SGESV()

    subroutine solve(A,bx)
        real(wp), dimension(:,:), intent(inout) :: A
        real(wp), dimension(:), intent(inout) :: bx
        integer N, nrhs, info,flag
        integer, dimension(:), allocatable :: ipiv
        
        N = size(A, dim=1)
        allocate(ipiv(N), stat= flag)
        if(flag /= 0 ) print *, "problem allocating ipiv"
        nrhs = 1
        call SGESV(N, nrhs, A, N, ipiv, bx, N, info)
        if(info /= 0) print *, "problem when using sgesf"
        deallocate(ipiv)

    end subroutine
    
end module solver