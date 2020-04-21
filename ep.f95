PROGRAM ep
    implicit none
    integer n,i,j
    real, allocatable :: v(:), w(:), D(:), L(:), x(:), y(:), z(:), b(:)
    n=5
    allocate(v(n),w(n-1),D(n),L(n-1), b(n))
     v = 0.; w = 0.; D = 0.; L = 0.; b = 1.

    !ARMAZENAMENTO DAS DIAGONALES
    v = (/4., 5., 3., 2., 3./); w=(/2., 1., 1., 1./)

    !FATORACAO DA MATRIZ
    call fatoracao_LDL(v,w,L,D,n)

    !SOLUCAO DO SISTEMA EQUIVALENTE
    x = sol_LDL(L,D,b,n)

    print*, x

    contains

    subroutine fatoracao_LDL(Adiag,Asub,L,D,n)
        implicit none
        integer :: i, n
        real :: Adiag(:), Asub(:), D(:), L(:) !Adiag DIAGONAL PRINCIPAL, Asub SUBDIAGONAL

        D(1) = Adiag(1)
        L(1)=Asub(1)/D(1)

        do j = 2, n-1
            D(j) = Adiag(j) - D(j-1)*L(j-1)**2
            L(j) = Asub(j)/D(j)
        end do

        D(n) = Adiag(n) - D(n-1)*L(n-1)**2
    end subroutine

    function sol_LDL(L,D,b,n)
        implicit none
        integer :: i,n
        real, allocatable :: x(:), y(:), z(:),sol_LDL(:)
        real :: L(:), D(:), b(:)
        allocate(x(n), y(n), z(n), sol_LDL(n))
        x = 0.; y = 0.; z = 0.; sol_LDL = 0.

        !RESOLVEMOS O SISTEMA Lz=b
        z(1) = b(1)
        do i = 1, n
            z(i) = b(i) - L(i-1)*z(i-1)
        end do

        !RESOLVEMOS O SISTEMA Dy=z
        do i = 1, n
            y(i) = z(i)/D(i)
        end do

        !RESOLVEMOS O SISTEMA L'x=y
        x(n) = y(n)
        do i = n-1, 1, -1
            x(i) = y(i) - L(i)*x(i+1)
        end do
        sol_LDL = x
    end function

END PROGRAM
