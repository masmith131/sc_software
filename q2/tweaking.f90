program tweaking 
    use siqrd_solver
    implicit none 

    integer :: num_args ! number of arguments passed in command line 
    character(len = 32) :: arg1, arg2, arg3, arg4 ! arguments of command line 
    real(wp) :: argT   ! simulation horizon
    real(wp) :: deltaB !delta_Beta value for FD approximation
    integer :: argN    ! N +1:  number of grid points in time interval [0,T]
    integer :: maxInf  !the target = Maximal number of Infections
    real(wp), dimension(7) :: input ! stores the input read from 'parameters.in'
    real(wp) optiB !optimal Beta
    integer iter 

    ! variables for timing 
    ! real(wp), dimension(9) :: times
    real(wp) :: t1, t2
    ! integer i 
    ! real(wp) time

    ! get N,T,deltaB, maxInf from command line
    num_args = command_argument_count()
    if(num_args /= 4 ) then 
        print *, "ERROR: wrong number of arguments provided"
        stop
    endif
    call get_command_argument(1,arg1)
    call get_command_argument(2,arg2)
    call get_command_argument(3,arg3)
    call get_command_argument(4,arg4)
    read(arg1,*)argN
    read(arg2,*)argT
    read(arg3,*)deltaB
    read(arg4,*)maxInf
   
    ! getting parameters and s0, i0 from the input file 
    open(unit = 32,file = 'parameters.in')
    read(32,*) input
    close(unit = 32)

    ! set parameters 
    call setting_parameters(input(:5), input(6:), argT, argN)

    ! WHEN DOING TIMINGS
    ! do i = 1, 9
    !     optiB = input(1)
    !     call cpu_time(t1)
    !     call newton(optiB)
    !     call cpu_time(t2)
    !     times(i) = t2-t1
    !     print *, optiB
    !     print *, iter 
    ! enddo
    ! print *, times
    ! time =  median(times)
    ! print *, time

    call cpu_time(t1)
    optiB = input(1)
    call newton(optiB)
    print *, "number of Newton iterations", iter
    print *, "optimal beta: ", optiB
    call cpu_time(t2)
    t2 = t2 -t1
    print *,"cpu_time: ", t2

    contains 

    ! WHEN DOING TIMINGS
    ! function median(times) result(m)
    !     real(wp), dimension(9) :: times
    !     integer i, j
    !     real(wp) temp, m
    !     ! Simple bubble sort to order the times array
    !     do i = 1, 8
    !         do j = i+1, 9
    !             if (times(i) > times(j)) then
    !                 temp = times(i)
    !                 times(i) = times(j)
    !                 times(j) = temp
    !             end if
    !         end do
    !     end do
    !     m = times(5)

    ! end function 

    !----------------------------------------------------------------------------------------------------------------
    ! f = fun(b)
    ! function of b that will be used for the newton method 
    ! function computes max Ik + Qk by calling iq_max
    ! using the method 'meth' 
    ! but can be changed to different function of beta 
    ! ----------------------------------------------------------------------------------------------------------------
    function fun(b) result(f)
        real(wp), intent(in) :: b  ! beta 
        real(wp) :: f  ! will contain the result F(beta)
        character :: meth = 'h' ! method used to solve siqrd system ('f': forward Euler, 'b': backward euler, 'h': heun)    
        ! call a function of beta here 
        f = iq_max(b, meth)
    end function 

    !---------------------------------------------------------------------------------------------------
    ! newton(B) performs Newton's method on variable beta 
    ! for function of beta: fun(b)
    ! result in B 
    !----------------------------------------------------------------------------------------------------
    subroutine newton(B)
        real(wp), intent(inout) :: B ! variable beta 

        real(wp) :: B_next ! Beta at next step of newton 
        real(wp) :: f_b ! F(beta)
        real(wp) :: f_bd ! F(beta + delta_beta)
        real(wp) :: dfdb ! finite difference approx of F'(beta)
        integer, parameter :: max_it = 200       ! maximum number of iterations to avoid infinite loops
        real(wp), parameter :: rtol = 1e-15_wp   ! relative tolerance for stopping criteria 
        real(wp), parameter :: eps = 1e-15_wp    ! shouldn't divide by smaller value 

        do iter = 0, max_it       
            f_b = fun(B)          ! evaluate function at beta 
            f_bd = fun(B + deltaB)    ! evaluate funcion at beta + delta_beta  
            dfdb = (f_bd - f_b) / deltaB ! get F'(beta) thanks to finite difference approximation 
            if(dfdb < eps) then 
                print *, "Derivative(denominator) is too small"
                exit 
            endif 
            B_next = B - (f_b - maxInf)/dfdb   !One step of Newton's Method 

            !stop criteria 
            if(abs(B_next - B) < rtol * abs(B)) then 
                B = B_next
                exit 
            endif 

            B = B_next 
        enddo 
    end subroutine 
    
end program