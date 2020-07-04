Program Rayleigh
!==================================================================
    !Apresentado por Deyanira Ribero Pineda
    !N° USP 11549649
    !Tarefa Computacional Splines
    !==================================================================
    implicit none
    integer :: n, item
    real(kind=8), parameter :: pi = 4.d0*Atan(1.d0)
    real(kind=8) :: lambda, h, a = 0.d0, b = 1.d0, erro
    real(kind=8), allocatable :: x(:), phi(:,:), phi2(:,:), diag(:), subdiag(:), d(:), sol(:), u(:,:)
    
    call Parametros(n,lambda,x,h,item)
    call spline(x, h, phi)
    call Amontagem(diag, subdiag, x, h, lambda)
    call dMontagem(d,x,h,lambda, item)
    call fatoracao_LDL(diag,subdiag)
    sol = sol_LDL(subdiag,diag,d)
    call Resultado(u,erro,item)
    call Exportar(u)
    deallocate(u)
    
    contains
    
    subroutine Parametros(n,lambda,x,h,item)
        implicit none
        
        integer :: n, item, i
        real(kind=8), allocatable :: x(:)
        real(kind=8) :: lambda, h
        
        print*, 'Digite o numero de nos interiores: '
        read*, n

        write(*,*) "Qual exercicio deseja rodar (1, 2 ou 3): "
        read*, item

        if (item .eq. 1) then
            print*, 'Digite lambda (deve ser positivo): '
            read*, lambda
        elseif (item .eq. 2) then
            lambda = 0.d0
        else
        !Aprensenta-se o lambda para o item 3 e para o teste 4 dependendo do caso que precise rodar 
            lambda = pi**2
            !lambda = 1.d0
        endif

        allocate(x(n+2), phi(2, 2),phi2(2, 2), diag(n), subdiag(n-1),d(n))
        h = (b - a) / (n + 1.d0); phi = 0.d0; x = 0.d0

        do i = 1, n+2
            x(i) = a + (i - 1.d0) * h
        enddo
        i = 0
    end subroutine
    
    subroutine Exportar(u)
        implicit none
        
        integer :: k, j
        real(kind=8) :: u(:,:)
        open(unit=10, File = "resultados.txt")
        k = 0; j = 0; k = size(u(:,1))
        
        do j = 1, k
            write(10,*)u(j,1), u(j,2), u(j,3)
        enddo
        
        close(10)
        
    end subroutine
    
    subroutine Resultado(u,E,item)
        implicit none
        
        integer :: item, k, j
        real(kind=8) :: E
        real(kind=8), allocatable :: u(:,:)
        
        k = 0; j = 0
        if ((item .eq. 2) .or. (item .eq. 3)) then
            k = 10*n
        allocate(u(k+1,3))
        
        do j = 1, k + 1
            u(j,1) = float((j-1)) / float(k)
            u(j,2) = exata(u(j,1),item)
            u(j,3) = aprox(u(j,1))
        enddo

        E = norm(u(:,2),u(:,3))       
        else
            allocate(u(n+2,3))
            u(:,1) = x
            do j = 1, n+2
                u(j,2) = exata(u(j,1),item)
                u(j,3) = aprox(u(j,1))
            enddo
            E = norm(u(:,2),u(:,3))
        endif
       ! print*,u(:,1)
       print*, E
    end subroutine
    
    
    !================================================
    !função chapéu
    !================================================
    subroutine spline(x, h, phi)
        implicit none
        integer :: ind
        real (kind=8) :: x(:), h, phi(:,:)

        phi = 1.d0; ind = 1
        phi(1,:) = (/-x(ind), 1.d0/)
        phi(2,:) = (/x(ind+2), -1.d0/)
        phi = (1.d0 / h) * phi

    end subroutine
    !========================================
        !Traslacao dos splines
    !========================================
    subroutine Translacao(p, q, r)
        implicit none
        real (kind=8) :: p(:,:), q(:,:), r

        q = p

        q(1,1) = -1.d0*(n + 1.d0) * r
        q(2,1) = 2.d0 + (n + 1.d0) * r

    end subroutine
    !=========================================
    !montagem da matriz do sistema normal
    !=========================================
    subroutine Amontagem(diag, subdiag, x, h, lambda)
        implicit none
        integer :: i,ind
        real(kind=8) :: diag(:), subdiag(:), h, lambda,x(:),aux
        aux = 1.d0 / h**2
        
        do i = 1,n
            ind = i + 1
            diag(i) = aux * (-x(ind-1) - (lambda / 3.d0)*(x(ind-1)**3) - lambda*x(ind-1)*(x(ind)**2) &
            + lambda*(x(ind-1)**2)*x(ind) + x(ind+1) - lambda*x(ind)*(x(ind+1)**2) + lambda*(x(ind)**2)*x(ind+1) &
            +(lambda / 3.d0)*x(ind+1)**3)
            ind = 0
        enddo
        i=0
        do i = 1,n - 1
            ind = i + 1
            subdiag(i) = aux * (x(ind+1) - x(ind+2) + (1.d0 / 6.d0) * lambda  * (-3.d0 * x(ind+1) * x(ind+2)**2 &
            + 3.d0 * (x(ind+1)**2) * x(ind+2) - x(ind+1)**3 + x(ind+2)**3))
            ind = 0
        enddo
    end subroutine

    !====================================
            !montagem do vetor d
    !====================================
    subroutine dMontagem(d,x,h,lambda,item)
    implicit none
        integer :: i,ind, item
        real(kind=8) :: d(:), h, lambda, x(:), aux
        aux = 1.d0 / h
        ind = 0; i = 0
        if (item .eq. 1) then
            do i = 1,n
                ind = i + 1
                d(i) = aux * ((-1.d0/ 24.d0) * lambda * x(ind-1)**4 - lambda * (1.d0/ 4.d0) * x(ind)**4 &
                - (1.d0 / 24.d0) * lambda * x(ind+1)**4 + (1.d0 / 12.d0) * lambda * x(ind-1)**3 &
                + (1.d0 / 3.d0) * lambda * x(ind)**3 + (1.d0 / 12.d0) * lambda * x(ind+1)**3 &
                + (1.d0 / 2.d0) * x(ind-1)**2 + x(ind)**2 + (1.d0 / 2.d0) * x(ind+1)**2 &
                + (1.d0 / 6.d0) * lambda * x(ind-1) * (x(ind)**3) - (1.d0 / 4.d0) * lambda * x(ind-1) * (x(ind)**2) &
                - x(ind-1)*x(ind) + (1.d0 / 6.d0) * lambda * (x(ind)**3) * x(ind+1) &
                - (1.d0 / 4.d0) * lambda * (x(ind)**2) * x(ind+1) - x(ind) * x(ind+1))

            enddo
            ind = 0; i = 0
        elseif (item .eq. 2) then
            do i = 1,n
                ind = i + 1
                d(i) = -2.d0 * aux * (0.5d0*(x(ind-1)**4) + 3.d0*(x(ind)**4) + 0.5d0*(x(ind+1)**4) - (x(ind-1)**3) &
                - 4.d0*(x(ind)**3) - (x(ind+1)**3) + 0.5d0*(x(ind-1)**2) + (x(ind)**2) + 0.5d0*(x(ind+1)**2) &
                - 2.d0*x(ind-1)*(x(ind)**3) + 3.d0*x(ind-1)*(x(ind)**2) - x(ind-1)*x(ind) - 2.d0*(x(ind)**3)*x(ind+1) &
                + 3.d0*(x(ind)**2)*x(ind+1) - x(ind)*x(ind+1))
            enddo
        ind = 0; i=0
        elseif (item .eq. 3) then
            do i = 1,n
                ind = i + 1
                d(i) = 2.d0 * aux * (pi * (x(ind-1) - 2.d0 * x(ind) + x(ind+1)) * cos(pi*x(ind)) - sin(pi*x(ind-1)) &
                + 2.d0 * sin(pi*x(ind)) - sin(pi*x(ind+1)))
                !O seguinte d e o usado para o teste 4 
                !d(i)=aux*(-6.d0*(x(ind)**4) - x(ind-1)**4 - x(ind+1)**4 + 8.d0*(x(ind)**3) + 2.d0*(x(ind-1)**3) &
                !+ 2.d0*(x(ind+1)**3) + 24.d0*(x(ind)**2) + 12.d0*(x(ind-1)**2) + 12.d0*(x(ind+1)**2) &
                !+ 4.d0*x(ind-1)*(x(ind)**3) - 6.d0*x(ind-1)*(x(ind)**2) - 24.d0*x(ind-1)*x(ind) &
                !+ 4.d0*(x(ind)**3)*x(ind+1) - 6.d0*(x(ind)**2)*x(ind+1) - 24.d0*x(ind)*x(ind+1))
            enddo
        endif
    end subroutine

    !=============================================
            !Resolucao do sistema normal
    !=============================================

    !=============================================
        !Fatoracao LDL Para o sistema normal
    !=============================================
    subroutine fatoracao_LDL(diag,subdiag)
        implicit none
        integer :: j
        real(kind=8), allocatable :: aux1(:), aux2(:)
        real(kind=8) :: diag(:), subdiag(:)
        allocate(aux1(size(diag)), aux2(size(subdiag)))

        j = 0
        aux1(1) = diag(1)
        aux2(1)=subdiag(1)/aux1(1)

        do j = 2, n-1
            aux1(j) = diag(j) - aux1(j-1)*aux2(j-1)**2
            aux2(j) = subdiag(j)/aux1(j)
        end do

        aux1(n) = diag(n) - aux1(n-1)*aux2(n-1)**2

        diag=0.d0; subdiag=0.d0
        diag=aux1; subdiag=aux2

    end subroutine

    !=================================
        !Solucao do sistema normal
    !=================================
    function sol_LDL(subdiag,diag,d)
            implicit none
            integer :: j
            real(kind=8), allocatable :: c(:), y(:), z(:),sol_LDL(:)
            real(kind=8) :: subdiag(:), diag(:), d(:)
            allocate(c(n), y(n), z(n), sol_LDL(n))
            c = 0.d0; y = 0.d0; z = 0.d0; sol_LDL = 0.d0; j = 0

            !RESOLVEMOS O SISTEMA Lz=d
            z(1) = d(1)
            do j = 1, n
                z(j) = d(j) - subdiag(j-1)*z(j-1)
            end do
            j=0
            !RESOLVEMOS O SISTEMA Dy=z
            do j = 1, n
                y(j) = z(j)/diag(j)
            end do

            j=0
            !RESOLVEMOS O SISTEMA L'c=y
            c(n) = y(n)
            do j = n-1, 1, -1
                c(j) = y(j) - subdiag(j)*c(j+1)
            end do
            sol_LDL = c
        end function
!===================================================
    !solucao aproximada da eq. diferencial
!===================================================
    function PolEval(p,r)
        implicit none

        real(kind=8) :: p(:), r, PolEval

        PolEval = p(1) + p(2)*r
    end function

    function aprox(w)
        implicit none

        integer :: j,ii
        real(kind=8) :: aprox, w, soma

        soma = 0.d0; j = 0

        do j = 1, n
            call Translacao(phi,phi2,h*(j-1))
            if (x(j) .le. w .and. w .le. x(j+2)) then
                if (x(j) .le. w .and. w .le. x(j+1)) then
                    soma = soma + sol(j) * PolEval(phi2(1,:),w)
                elseif (x(j+1) .le. w .and. w .le. x(j+2)) then
                    soma = soma + sol(j) * PolEval(phi2(2,:),w)
                endif
            endif
            phi2 = 0.d0
        enddo
        aprox = soma
    end function
 !======================================================
    !SOLUCAO EXATA DA EQ. DIFERENCIAL
 !======================================================
    function exata(w, item)
        implicit none
        integer :: item
        real(kind=8) :: w, exata

        if (item .eq. 1) then
            exata = 0.5d0 * w * (1.d0 - w)
        elseif (item .eq. 2) then
            exata = (w**2) * (w - 1.d0)**2
        elseif (item .eq. 3) then
            exata = sin(pi*w)
            !A seguinte funcao exata e o caso do teste 4 no caso de que quer ou precisar rodar esse teste
           !  exata = 12.d0 * w * (1.d0 - w)
        endif
     end function
!========================================================
        !NORMA INFINITO
!========================================================
    function norm(z,v)
        implicit none
        integer :: j
        real(kind=8) :: z(:), v(:), norm, aux1, aux2
        aux1 = abs(z(1)-v(1))
        do j = 2, size(z)
            aux2 = abs(z(j)-v(j))
            norm = max(aux1,aux2)
            aux1 = norm
        end do
    end function
end program
