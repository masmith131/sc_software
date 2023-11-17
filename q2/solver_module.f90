module solver
    !use lapack
    implicit none

    integer, parameter :: s_p = kind(1.0) ! single precision 
    integer, parameter :: d_p = kind(1.0d0)  ! double precision 
    ! following line to change for different precision 
    integer, parameter, public :: wp=d_p ! set select_p to the desired precision (sp or dp)
    
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

        if(info /= 0) print *, "problem when using sgesf/dgesf", info
        deallocate(ipiv)

    end subroutine


    subroutine get_eigenvalues(A, ev)
        real(wp), dimension(5,5), intent(in) :: A !5x5 matrix for which we are computing the eigenvalues
        real(wp), dimension(5), intent(out) :: ev !array of real parts of eigenvalues

        integer :: info, lwork  !flag and LAPACK routine required parameter
        real(wp), dimension(5) :: wr, wi ! ! real and imaginary parts of eigenvalues
        real(wp), dimension(5*5) :: work
        lwork = 15
        call DGEEV('N', 'N', 5, A, 5, wr, wi, work, 1, work, 1, work, lwork, info)
        if (info /= 0) print *,"problem when using DGEEV", info
        print *, wi 
        print *, wr 
        ev = wr
    end subroutine 

    ! subroutine get_eigenvalues(A, ev)
    !     real(wp), intent(in), dimension(5,5) :: A
    !     real(wp), intent(out), dimension(5) :: ev

    !     real(wp), dimension(5*5) :: work
    !     real(wp), dimension(5) :: wr, wi  
    !     integer :: info, lwork
    !     lwork = 15
    
    !     call DGEEV('N', 'N', 5, A, 5, wr, wi, work, 1, work, 1, work, lwork, info)
    
    !     ! Check for success
    !     if (info /= 0) then
    !         print *, 'An error occurred: ', info
    !         stop
    !     end if
    
    !     ! Since A is real, we are only interested in the real part of the eigenvalues.
    !     ! Copy the real part of the eigenvalues to the output variable ev.
    !     ev = wr
    
    ! end subroutine get_eigenvalues
    
    
end module solver