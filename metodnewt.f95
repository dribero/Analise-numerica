program metod
    implicit none
    integer :: i
    !v e a diferenca do ponto p0 e pn
    !p0 e o chute inicial e pn e a nova iteracao

    real(kind=8) ::  P0, tol, iter, pn, v
    P0 = 1.d0; tol = 1e-5; iter = 100; i = 1

    do while (i <= iter)
        pn = P0 - f(P0) / der(P0)
        v = pn - P0

        if (abs(v) .lt. tol) then
            exit
        end if

        i = i + 1
        P0 = pn

    end do

    print*, pn, f(pn)

    contains

    function f(x)
        implicit none
        real(kind=8) :: x, f
        f = x**2 - 2
    end function

    function der(x)
        implicit none
        real(kind=8):: x, der
        der = 2*x
    end function


end program
