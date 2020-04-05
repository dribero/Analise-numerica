PROGRAM secantes
    implicit none
    integer :: i
    !P0 e P1 sao os chutes iniciais, TOLerancia, iter e o numero maximo de iteracoes, pn e a nova iteracao
    real(kind=16) ::  P0, P1, tol, iter, pn, v, exata, xold, xnew
    P1 = 1.d0; P0 = 2.d0; TOL = 1e-4; iter = 100; i = 1; exata=sqrt(2.d0); xold=P0

    do while (i <= iter)
        pn = P1 - (f(P1)*(P1-P0)) / (f(P1)-f(P0))
        v = pn - P1

        if (abs(v) .lt. tol) then
            exit
        end if

        i = i + 1
        P0 = P1
        P1 = pn

        xold=xnew
        xnew=pn
    end do

    print*, i, Pn, o(xold,xnew)

    contains

    function f(x)
        implicit none
        real(kind=16) :: x, f
        f = x**2 - 2
    end function

    function o(xold,xnew)
    implicit none
    real*16 :: xold, xnew, o, num, den
    num = abs(xnew - exata); den = abs(xold - exata)
    o = num / den
    end function

END PROGRAM
