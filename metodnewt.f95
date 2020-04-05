program metod
    implicit none
    integer :: i
    !v e a diferenca do ponto p0 e pn
    !p0 e o chute inicial e pn e a nova iteracao

    real(kind=8) ::  P0, tol, iter, pn, v, exata, xold, xnew
    P0 = 1.d0; tol = 1e-5; iter = 100; i = 1; exata = sqrt(2.d0); xold = P0

    do while (i <= iter)
        pn = P0 - f(P0) / der(P0)
        v = pn - P0

        if (abs(v) .lt. tol) then
            exit
        end if

        i = i + 1
        P0 = pn

        xold=xnew
        xnew=Pn
    end do

    print*, i, Pn, o(xold,xnew)

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

    function o(xold,xnew)
    implicit none
    real*8 :: xold, xnew, o, num, den
    num = abs(xnew - exata); den = abs(xold - exata)**2
    o = num / den
    end function

end program
