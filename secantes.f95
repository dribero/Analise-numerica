PROGRAM secantes
    implicit none
    integer :: i
    !P0 e P1 sao os chutes iniciais, tol e a tolerancia, iter e o numero maximo de iteracoes, pn e a nova iteracao
    real(kind=8) ::  P0, P1, tol, iter, pn, v
    P1 = 1.d0; P0 = 2.d0; tol = 1e-5; iter = 100; i = 1

    do while (i <= iter)
        pn = P1 - (f(P1)*(P1-P0)) / (f(P1)-f(P0))
        v = pn - P1

        if (abs(v) .lt. tol) then
            exit
        end if

        i = i + 1
        P0 = P1
        P1 = pn

    end do

    print*, pn, f(pn)

    contains

    function f(x)
        implicit none
        real(kind=8) :: x, f
        f = x**2 - 2
    end function

END PROGRAM
