module solver
    implicit none

    integer, parameter :: s_p = kind(1.0) ! single precision 
    integer, parameter :: d_p = kind(1.0d0)  ! double precision 
    ! following line to change for different precision 
    integer, parameter, public :: wp=d_p ! set select_p to the desired precision (sp or dp)
    integer, parameter :: s_ip = selected_int_kind(9)
    integer, parameter :: d_ip = selected_int_kind(15)
    integer, parameter, public :: wip = d_ip 
    
    contains

    !---------------------------------------------------------------------------------------------
    ! solve(A,bx) solves the linear system Ax = b
    ! bx contains b at function call 
    ! bx contains x after function was exectuted
    ! it uses the Lapack subroutine SGESV()
    !---------------------------------------------------------------------------------------------

    subroutine solve(A,bx)
        real(wp), dimension(:,:), intent(inout) :: A
        real(wp), dimension(:), intent(inout) :: bx
        integer N, nrhs, info,flag
        integer, dimension(:), allocatable :: ipiv
        
        N = size(A, dim=1)
        allocate(ipiv(N), stat= flag)
        if(flag /= 0 ) print *, "problem allocating ipiv"
        nrhs = 1

        if(wp == s_p) then
            call SGESV(N, nrhs, A, N, ipiv, bx, N, info)
        elseif(wp == d_p) then 
            call DGESV(N, nrhs, A, N, ipiv, bx, N, info) 
        else 
            print *, "Problem reconginzing the kind"
        endif


        if(info /= 0) print *, "problem when using sgesf"
        deallocate(ipiv)

    end subroutine
    
end module solver