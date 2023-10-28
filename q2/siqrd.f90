program siqrd
    use siqrd_solver 
    implicit none 

    ! following character contains the name of the method we will be using 
    character, parameter :: m = 'b' ! f = forward euler, b = backward euler, h= heun

    real, parameter :: step = T/N !time step between grid points
    real, dimension(N+1) :: grid !array of size N+1 containg grid points on [0,T]
    real, dimension(7) :: input ! stores the input read from 'parameters.in'
    real, dimension(5,N+1) :: sol ! stores the solution at each time step 
    integer i,j

    ! fill in array containg grid points 
    grid(1) = 0.0
    do i = 2, N+1
        grid(i) = grid(i-1) + step
    enddo

    ! getting parameters and s0, i0 from the input file 
    open(unit = 32,file = 'parameters.in')
    read(32,*) input
    close(unit = 32)

    ! xk = (s(t), i(t), q(t), r(t), d(t)), corresponds to one row of sol 
    ! filling in x0
    sol(:,1) = 0.0
    sol(1,1) = input(6) ! s0
    sol(2,1) = input(7) ! i0

    ! setting parameters for the numerical method method 
    call setting_parameters(input(:5))

    ! calling relevant method for every time step 
    do i = 2, N+1
        if(m == 'f') then 
            call forward(sol(:,i-1),sol(:,i))
        elseif(m == 'h') then 
            call heun(sol(:,i-1), sol(:,i))
        elseif(m == 'b') then 
            print *, "iteration",i
            call backward(sol(:,i-1),sol(:,i))
        else 
            print *, "method m is not recognized"
            exit
        endif
    enddo

    ! printing the solution vector for every time step 
    do j = 1, N+1
        write(*,'(es12.5, 1x, a)',advance='no') grid(j), " "
        do i = 1, 5
            write(*, '(es12.5, 1x, a)', advance='no') sol(i,j), " "
        enddo
        print *
    enddo


end program 