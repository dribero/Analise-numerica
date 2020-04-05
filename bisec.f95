program bisec
    implicit none
    integer :: iter, i
    real*16 :: a,b,tol,raiz,p,fa,fp,exata,xold,xnew ![a,b], TOLerancia, raiz

    a = 1.d0; b = 2.d0; tol = 1e-14; iter = 100; i = 0; p = (a + b) / 2.d0; exata = sqrt(2.d0); xold = p

    do while ((i .le. iter) .and. (abs(f(p)) .gt. tol) .and. ((b-a) / 2.d0 .gt. tol))

        i = i + 1

        if (f(a) * f(p) .gt. 0) then
            a = p
        else
            b = p
        endif

        p = (a + b) / 2.d0

        xold = xnew
        xnew = p

    enddo
print*, i,p, o(xold,xnew)

contains

function f(x)
    implicit none
    real*16 :: x, f
    f = (x**2) - 2.d0
end function

function o(xold,xnew)
    implicit none
    real*16 :: xold, xnew, o, num, den
    num = abs(xnew - exata); den = abs(xold - exata)
    o = num / den
end function
end program
