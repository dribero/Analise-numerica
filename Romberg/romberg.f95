program romberg
!==================================================================
    !Apresentado por Deyanira Ribero Pineda
    !N° USP 11549649
    !Tarefa Computacional Romberg
    !==================================================================
    implicit none
    real*8 :: a, b, eps = 1e-6
    real*8, allocatable :: R(:)
    character :: cont
    integer :: n = 6, teste, itmax = 100

    allocate(R(n+1))

    R = 0.d0
    call main()

    do
        print*, "Deseja continuar com outro teste? (s: sim, n: nao) "
        read*, cont
        if (cont == "s") then
            R = 0.d0
            call main()
        else
            exit
        endif
    enddo

    contains

    subroutine main()
        implicit none

        write(*,*) "Digite o teste a fazer (1: primeira integral, 2: segunda integral, 3: funcao de Runge)"
        read*, teste

        if (teste == 1) then
            a = 0.d0; b = 1.d0
        elseif (teste == 2) then
            a = 0.d0; b = 0.995d0
        else
            a = -5.d0; b = 5.d0
        endif

        call romb()
    end subroutine

    function f(x)
        implicit none
        real*8 :: f, x

        if (teste == 1) then
            f = x**2
        elseif (teste == 2) then
            f = 1.d0 / (1.d0 - x)
        else
            f = 1.d0 / (1.d0 + x**2)
        endif
    end function

    function trapz(t, i)
        implicit none
        real*8 :: trapz, t, h
        integer :: i, j, k

        trapz = 0.d0; h = (b - a) / 2**i

        if (i == 0) then
            trapz = (b - a) / 2.d0
            trapz = trapz * (f(a) + f(b))
        else
            k = 2**(i-1)
            do j = 1, k
                trapz = trapz + f(a + (2*j - 1)*h)
            enddo
            trapz = (h*trapz) + t / 2.d0
        endif
    end function

    subroutine romb()
        implicit none
        real*8 :: erro, aux
        integer :: i, k, j

        i = 0; aux = 0.d0

        do
            k = 1; i = i + 1
            R(1) = trapz(aux,(i-1))
            100 format (1X, A7, 1X, I2, 1X, A2)
            110 format (1X, A50, A50, A50, A30)
            write(*,100) "Tabela ", i, ": "
            write(*,110) "--------------------------------------------------","--------------------------------------------------" &
            ,"--------------------------------------------------","------------------------------"
            print*, R(1)
            R(2:) = 0.d0

            if (k .eq. 1) then
                    aux = R(k)
                endif

            do k = 2, n+1
                R(k) = trapz(R(k-1),k+i-2)

                do j = 1, k
                    R(k-j) = (R(k-j+1) - R(k-j)) / ((4**j) - 1.d0)
                    R(k-j) = R(k-j) + R(k-j+1)
                enddo

                print*, R(k:1:-1)
            enddo

            if (abs(R(n) - R(n+1)) .le. eps*abs(R(n+1)) .or. i .ge. itmax) then
                print*
                print*, "O erro relativo e: ", abs(R(n) - R(n+1))
                print*
                if (teste == 1) then
                        print*, "O erro absoluto é: ", abs(R(n+1)-(1.d0/3.d0))
                    elseif (teste == 2) then
                        print*, "O erro absoluto é: ", abs(R(n+1)-log(200.d0))
                    else
                        print*, "O erro absoluto é: ", abs(R(n+1)-2.d0*atan(5.d0))
                endif

                exit
            endif
            print*
        enddo

    end subroutine

end program
