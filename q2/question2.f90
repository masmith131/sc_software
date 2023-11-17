program siqrd
    use siqrd_solver 
    implicit none 

    character, parameter :: m = 'h' !name of the method :  f = forward euler, b = backward euler, h= heun
    real(wp) :: step !time step between grid points
    real(wp) :: argT! simulation horizon
    integer :: argN ! N +1:  number of grid points in time interval [0,T]
    integer :: num_args ! number of arguments passed in command line 
    character(len = 32) :: arg1, arg2

    real(wp), allocatable :: sol(:,:) ! will contain solution for all time steps
    real(wp), allocatable :: grid(:)  ! will contain all time steps

    real(wp), dimension(7) :: input ! stores the input read from 'parameters.in'
    integer i,j
    integer flag1, flag2

    ! getting T and N from the command line
    num_args = command_argument_count()
    if(num_args /= 2 ) then 
        print *, "ERROR: wrong number of arguments provided"
        stop
    endif
    call get_command_argument(1,arg1)
    call get_command_argument(2,arg2)
    read(arg1,*)argN
    read(arg2,*)argT

    allocate(grid(argN+1),stat = flag1) !array of size N+1 containg grid points on [0,T]
    allocate(sol(5, argN+1), stat = flag2) ! stores the solution at each time step 
    if(flag1 /= 0 .or. flag2 /= 0) print *, "Problem with allocation of grid or sol"
    step = argT/real(argN, wp)

    ! fill in array containg grid points 
    grid(1) = 0.0
    do i = 2, argN+1
        grid(i) = grid(i-1) + step
    enddo

    ! getting parameters and s0, i0 from the input file 
    open(unit = 32,file = 'parameters.in')
    read(32,*) input
    close(unit = 32)

    ! setting parameters for the numerical method method 
    call setting_parameters(input(:5), input(6:), argT, argN)

    ! xk = (s(t), i(t), q(t), r(t), d(t)), corresponds to one row of sol 
    ! filling in x0
    sol(:,1) = 0.0_wp
    sol(1,1) = S0 ! s0
    sol(2,1) = I0 ! i0

    ! calling relevant method for every time step 
    do i = 2, argN+1
        if(m == 'f') then 
            call forward(sol(:,i-1),sol(:,i))
        elseif(m == 'h') then 
            call heun(sol(:,i-1), sol(:,i))
        elseif(m == 'b') then 
            call backward(sol(:,i-1),sol(:,i))
        else 
            print *, "method m is not recognized"
            exit
        endif
    enddo

    ! printing the solution vector for every time step 
    do j = 1, argN+1
        write(*,'(es12.5, 1x, a)',advance='no') grid(j), " "
        do i = 1, 5
            write(*, '(es12.5, 1x, a)', advance='no') sol(i,j), " "
        enddo
        print *
    enddo

    deallocate(sol)
    deallocate(grid)


end program 