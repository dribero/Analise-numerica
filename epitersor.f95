Program EPITERSOR
    !==================================================================
    !Apresentado por Deyanira Ribero Pineda
    !N° USP 11549649
    !Tarefa Computacional Métodos iterativos para Sistemas Lineares
    !==================================================================
    implicit none
	integer, parameter :: n = 50
	integer :: dim, dim1, dim2, ind
	integer :: i, j, maxiter = 5000, K = 100, chute
	integer, allocatable :: iter(:)
	real(kind=8) :: wi, tol = 1e-4
	real(kind=8), dimension(n) :: borda = 0.d0
	real(kind=8), allocatable :: subdiag(:), b(:), xinic(:), x(:,:), aux(:)

	dim = (n - 1)**2; dim1 = n*(n - 2)

	do j = 1, n
		borda(j) = ((6.d0 * j) / n) - 3.d0
	enddo

	allocate(subdiag(dim1), b(dim), xinic(dim), x(K + 1,dim), aux(dim), iter(K + 1))

    !Escolha da aproximacao inicial

	write(*,*)"Escolha o tipo de chute inicial (1: Nulo ou 2: Numeros aleatorios entre 0 e 1) "
	read*, chute

	if (chute .eq. 1) then
		xinic = 0.d0
	else
		do i = 1, dim
			xinic(i) = rand()
		enddo
	endif

	call montagem(subdiag, b)

    !arquivo de texto que armazena a solucao com o omega dado pelo menor numero de iteracoes
	Open(Unit = 130, File = "solucao.txt", ACCESS = "SEQUENTIAL")


	do i = 1, K + 1
		wi = 1.d0 + (real(i - 1) / real(K))
		x(i,:) = xinic
		call sor(x(i,:), subdiag, b, aux, wi, iter(i))
		write(*,100) i, 1.d0 + (real(i - 1) / real(K)), iter(i)
		100 format (1X,"Iteracao terminada com w(", I4, 1X, ") = ", F8.5, 1X, ", Iter. = ", I5)
	enddo

    !ind é o indice correspondente ao wi que produz menor numero de iteracoes
	ind = miniter(iter)

	print*
	write(*,110) ind, 1.d0 + (real(ind - 1) / real(K))
	110 format (1X,"w(", I2, 1X, ") = ", F8.5, 1X)
	print*
	write(*,120) iter(ind)
	120 format (1X,"Numero de iteracoes = ", I5)
	print*


    ! A solucao no arquivo de texto mencionado anteiormente
	do i = 1, dim
		write(130,*) x((ind),i)
	enddo

    close(Unit = 130)


	contains

    !a norma infinito da diferença de dois vetores
	function norm(u,v)
		implicit none
		integer :: i = 0
        real(kind=8) :: u(:), v(:), norm, aux1, aux2

        aux1 = abs(u(1)-v(1))
        do i = 2, size(u)
            aux2 = abs(u(i)-v(i))
            norm = max(aux1,aux2)
            aux1 = norm
        end do
	end function

    !Na seguinte funcao se calcula o indice do menor numero de iteracoes
	function miniter(iter)
		implicit none
		real(kind=8) :: aux1, aux2
		integer :: iter(:), i = 0, j = 1, miniter

		aux1 = iter(1)

		do i = 2, K + 1
			aux2 = iter(i)

			if (aux1 .gt. aux2) then
				aux1 = aux2
				j = i
			endif
		enddo

		miniter = j

	end function

        !Nesta subroutine se mostra como se armazena a informacao da matriz,
        ! e a informacao do vetor b.
        !Levando em conta que nao foi necessario armazenar a diag principal
        !somente foi armazenada uma subdiagonal
	subroutine montagem(subdiag, b)
		implicit none
		integer :: i = 0, l = 0
		real(kind=8) :: subdiag(:), b(:)

		subdiag = -1.d0; b = 0.d0

		do i = 1, n - 2
			l = (n - 1) * i
			subdiag(l) = 0.d0
		enddo

		b(1) = -6.d0; b(n - 1) = -3.d0 + borda(1); b(dim) = 2*borda(n - 1); b(dim - n + 2) = -3.d0 + borda(1)
		b(2:n - 2) = -3.d0; b(dim - n + 3:dim - 1) = borda(2:n - 2)

		do i = 1, n - 3
			l = n - 1
			b(i*l + 1) = -3.d0
			b(l*(i + 1)) = borda(i + 1)
		enddo

	end subroutine

        !O metodo SOR adaptado ao sistema do problema
	subroutine sor(x, subdiag, b, aux, w, cont)
		implicit none
		integer :: i = 0, j = 0, l, cont
		logical :: bool
		real(kind=8) :: subdiag(:), b(:), x(:), aux(:), w

		aux = 0.d0; l = 0; i = 0; bool = .true.

		do while ((l .lt. maxiter) .and. (bool .eqv. .true.))

			l = l + 1; aux = x; cont = l

			x(1) = ((1.d0 - w) * x(1)) + ((w / 4.d0) * (b(1) - subdiag(1)*x(2) + x(n)))

			do i = 2, n - 1
				x(i) = ((1.d0 - w) * x(i)) + ((w / 4.d0) * (b(i) - subdiag(i - 1)*x(i - 1) - subdiag(i)*x(i + 1) + x(i + n - 1)))
			enddo

			j = 0

			do i = n, dim - n + 1
				j = j + 1
				x(i) = ((1.d0 - w) * x(i)) + ((w / 4.d0) * (b(i) + x(j) + x(i + n - 1) - subdiag(i - 1)*x(i - 1) - subdiag(i)*x(i + 1)))
			enddo

			do i = dim - n + 2, dim - 1
				j = j + 1
				x(i) = ((1.d0 - w) * x(i)) + ((w / 4.d0) * (b(i) + x(j) - subdiag(i - 1)*x(i - 1) - subdiag(i)*x(i + 1)))
			enddo

			x(dim) = ((1.d0 - w) * x(dim)) + ((w / 4.d0) * (b(dim) - subdiag(dim1)*x(dim - 1) + x(dim - n + 1)))

			if (norm(aux,x) .le. tol) then
				bool = .false.
			endif

		enddo

	end subroutine



end Program
