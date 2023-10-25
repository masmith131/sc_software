program siqrd

    use models 

    ! following string contains the name of the method we will be using 
    ! f = forward euler, b = backward euler, h= heun
    character, parameter :: m = 'f'
    ! N +1:  number of grid points in time interval [0,T]
    integer, parameter:: N = 150
    ! T: simulation horizon 
    real, parameter :: T = 30
    ! input: will contain what is read from 'parameters.in'
    ! xk = (s(t), i(t), q(t), r(t), d(t))
    ! grid: array of size N+1 containg grid points on [0,T]
    ! step: time step between grid points
    real :: input(7), grid(N+1), step, sol(N+1,5)
    integer i,j 

    ! getting parameters and s0, i0 from the input file 
    open(unit = 32,file = 'parameters.in')
    read(32,*) input
    close(unit = 32)

    sol(1,:) = 0.0
    sol(1,1) = input(6)
    sol(1,2) = input(7)

    ! fill in array containg grid points 
    step = T/N
    grid(1) = 0.0
    do i = 2, N+1
        grid(i) = grid(i-1) + step
    enddo

    
    do i = 2,N+1
        if(m == 'f') then 
            call forward(input(:5), sol(i-1,:), T, N, sol(i,:))

        elseif(m == 'h') then 
            call heun(input(:5),sol(i-1,:), T, N, sol(i,:))
        elseif(m == 'b') then 
            print *, "BE not yet implemented ! "
            exit
        else 
            print *, "method m is not recognized"
            exit
        endif
    enddo


    do i = 1, N+1
        do j = 1, 5
            write(*, '(es12.5, 1x, a)', advance='no') sol(i,j), " "
        enddo
        print *
    enddo


end program 