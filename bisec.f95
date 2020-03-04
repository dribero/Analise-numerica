program bisec
    implicit none
    integer :: iter, i
    real*8 :: a,b,tol,raiz,p,fa,fp ![a,b], TOLerancia, raiz


    a = 1.d0;b = 2.d0; tol = 1e-5; iter = 100; i = 1

    fa = f(a)

    do while (i .le. iter)
        p = a + (b - a) / 2.d0
        fp = f(p)

        if ((fp .eq. 1e-5) .or. (((b - a) / 2) .lt. tol)) then
            exit
        endif

        i = i + 1

        if (fa * fp .gt. 1e-5) then
            a = p
            fa = fp
        else
            b = p
        endif
    enddo

    print*, p

contains

function f(x)
    implicit none
    real*8 :: x, f
    f = (x**2) - 2
end function

end program
