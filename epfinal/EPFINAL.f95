program ProjetoFinal

!==================================================================
    !Apresentado por Deyanira Ribero Pineda
    !N° USP 11549649
    !EP FINAL
!==================================================================
!NP é o número de pontos da malha para fazer o gráfico da função resultante, bem como o cálculo da norma
!Inf é a condição de contorno para x = 0
!Sup é a condição de contorno para x = 1
!no caso de colocar uma condicao diferente de contorno nao homogeneo,
!você deve definir a derivada da função k(x) na função "Dk" definida neste código
    implicit none
    integer, parameter :: single = selected_real_kind(p=4)
    integer, parameter :: double = selected_real_kind(p=8)
    integer, parameter ::  NP=10000
    real(kind=double) ::  Inf , Sup, h
    real(kind=double), allocatable, dimension(:) :: y, x1, y1, x, b
    real(kind=double), allocatable, dimension(:,:) :: A
    integer :: i, j, linha, n, m, CondHom
    print*, "Escreva o numero de nos interiores (n = 15, 31, 63, 127 ou 255): "
    read*, n

    print*, "Qual tipo de splines deseja utilizar (1: lineares ou 2: cubicos)"
    read*, m

    m = 2*m
    print*, "Deseja condicoes de fornteira homogeneas? (1: sim, 0: nao): "
     read*, CondHom

     if (CondHom .eq. 1) then
        Inf = 0._double; Sup = 1._double
     else
        print*, "Escreva o valor do extremo inferior Inf: "
        read*, Inf
        print*
        print*, "Escreva o valor do extremo superior Sup: "
        read*, Sup
     endif

    allocate ( x(n+m-2), b(n+m-2))

    if (m.EQ.2) then
        allocate(A(3,n))
        do i = 1,3
            do j = 1, n
                A(i,j)=0._double
            end do
        end do
    elseif (m.EQ.4) then
        allocate(A(7,n+2))
        do i = 1,7
            do j = 1, n+2
                A(i,j)=0._double
            end do
        end do
    end if

    h=1._double/(n+1._double)

    if (m.EQ.2) then

            do i = 1, n
                linha = 2
                do j =i, min(i+1,n)
                    A(linha,i)=ProdutInterL(j,i)
                    linha=linha-1
                end do
                b(i)=ProdutInterF(i)
            end do
            A(3,:)=A(1,:)

        call Gauss(A,b,x,n,3)

        allocate(x1(NP+1))
        allocate(y(NP+1))
        allocate(y1(NP+1))

        h=1._double/(NP+1._double)
        x1(1)=0._double

        do i = 2, NP+1
            x1(i)=x1(i-1)+h
        end do

        do i = 1, NP+1
           y(i)=0._double
        end do

        y1=y

        h=1._double/(n+1._double)
        !================================================================================
        !Avaliacao da exata e aprox
        !===============================================================================
            do i = 1, NP+1
                y(i)=u(x1(i))
                do j = 1, n
                    y1(i)=y1(i)+x(j)*BSplinLin(x1(i),j)
                end do
                if ((Inf.NE.0._double).OR.(Sup.NE.1._double)) then
                    y1(i) = y1(i) + Inf + (Sup-Inf)*x1(i)
                end if
            end do

            call OUTPUTS(x1,y,NP+1,"exata.dat")
            call OUTPUTS(x1,y1,NP+1,"aprox.dat")
        !================================================================================
        !print do erro maximo
        !================================================================================
            print *, maxval(y-y1)


    else

            do i = 1, n + 2
                linha = 4
                do j =i, min(i+3,n+2)
                    A(linha,i)=ProdutInterLSC(j,i)
                    linha=linha-1
                end do
                b(i)=ProdutInterFSC(i)
            end do
            A(5,:)=A(3,:)
            A(6,:)=A(2,:)
            A(7,:)=A(1,:)
            print*,A(4,:)

        call Gauss(A,b,x,n+2,7)

        allocate(x1(NP+1))
        allocate(y(NP+1))
        allocate(y1(NP+1))

        h=1._double/(NP+1._double)
        x1(1)=0._double

        do i = 2, NP+1
            x1(i)=x1(i-1)+h
        end do

        do i = 1, NP+1
           y(i)=0._double
        end do

        y1=y

        h=1._double/(n+1._double)
        !===============================================================================
        !Avaliamos exata e Aprox
        !===============================================================================
        do i = 1, NP+1
            y(i)=u(x1(i))
            do j = 1, n + 2
                y1(i)=y1(i)+x(j)*BSC(x1(i),j)
            end do
            if ((Inf.NE.0._double).OR.(Sup.NE.1._double)) then
                y1(i) = y1(i) + Inf + (Sup-Inf)*x1(i)
            end if
        end do

            call OUTPUTS(x1,y,NP+1,"exata.dat")
            call OUTPUTS(x1,y1,NP+1,"aprox.dat")

        !============================
        !print do erro maximo
        !============================
            print *, maxval(y-y1)

    end if


    contains

        function f(x)
        real(kind=double) :: f, x
        if ((Inf.EQ.0._double).AND.(Sup.EQ.1._double)) then
          !f = 12._double*x*(1._double-x)-2!Test1
         f = x + (2._double - x)*exp(x) !Test2
         !f = (2.d0-(x+1.d0)*(x+1.d0))*exp(1.d0)*log(2.d0)-2.d0*exp(x) !Test3

        else
         !f = -12._double*x**2 + 12._double*x - 2 !Test1
         f = exp(x)+1._double+exp(x)*((x-1._double)*(exp(-x)-1._double) + Inf + (Sup - Inf)*x)!Test2
         !f = (2._double-(x+1._double)*(x+1._double))*exp(1._double)*log(2._double)-2._double*exp(x) !Test3

        end if

    return
    end function f

     function q(x)
        real(kind=double) :: q, x
        !q = 0._double !Test1
        q = exp(x) !Test2
        !q = x+2.d0 !Test3
        return
    end function

        function k(x)
        real(kind=double) :: k, x
        !k = 1._double !Test1
        k = exp(x) !Test2
        !k = x+1.d0 !Test3
        return
    end function k

        function Dk(x)
        real(kind=double) :: Dk, x
        !Dk = 0._double !Test1
        Dk = exp(x) !Test2
        !Dk = 1.d0 !Test3
        return
    end function Dk


        function u(x)
        real(kind=double) :: u, x
        if ((Inf.EQ.0._double).AND.(Sup.EQ.1._double)) then
            !u = (x**2)*(x-1)**2 !Test1
            u = (x-1._double)*(exp(-x)-1._double) !Test2
            !u = exp(x)*log(x+1.d0)-(exp(1.d0)*log(2.d0))*x !Test3

        else
            !u = (x**2)*(x-1)**2 + Inf + (Sup - Inf)*x !Test1
            u = (x-1._double)*(exp(-x)-1._double) + Inf + (Sup - Inf)*x !Test2
            !u = exp(x)*log(x+1.d0)-(exp(1.d0)*log(2.d0))*x !Test3
        end if

    end function

     function FProdIntSplinLin(x,Bi,Bj)

        real(kind=double) :: x, FProdIntSplinLin
        integer :: Bi, Bj
        FProdIntSplinLin=k(x)*DBSplinLin(x,Bi)*DBSplinLin(x,Bj)
        FProdIntSplinLin=FProdIntSplinLin + q(x)*BSplinLin(x,Bi)*BSplinLin(x,Bj)
        return
    end function FProdIntSplinLin

     function FProdIntLSC(x,Bi,Bj)
        real(kind=double) :: x, FProdIntLSC
        integer :: Bi, Bj
        FProdIntLSC=k(x)*BSCub1(x,Bi,2)*BSCub1(x,Bj,2)
        FProdIntLSC=FProdIntLSC + q(x)*BSCub1(x,Bi,1)*BSCub1(x,Bj,1)
        return
    end function FProdIntLSC



    function BSplinLin(x,Bn)

        implicit none
        real(kind=double) :: x, BSplinLin
        integer :: Bn
        BSplinLin = x - (Bn-1)*h
        if ((BSplinLin.LT.0._double).OR.(BSplinLin.GT.2._double*h)) then
            BSplinLin=0._double
            return
        end if
        if ((BSplinLin.GE.0._double).AND.(BSplinLin.LE.h)) then
            BSplinLin = BSplinLin/h
            return
        else
            BSplinLin = 2._double - BSplinLin/h
            return
        end if
    end function BSplinLin

    function DBSplinLin(x,Bn)

        implicit none
        real(kind=double) :: x, DBSplinLin
        integer :: Bn
        DBSplinLin = x - (Bn-1)*h
        if ((DBSplinLin.LE.0._double).OR.(DBSplinLin.GE.2._double*h)) then
            DBSplinLin=0._double
            return
        end if
        if ((DBSplinLin.GT.0._double).AND.(DBSplinLin.LT.h))then
            DBSplinLin = 1._double/h
            return
        elseif ((DBSplinLin.GT.h).AND.(DBSplinLin.LT.2._double*h)) then
            DBSplinLin = -1._double/h
            return
        end if
    end function DBSplinLin


    function BSCub1(x,Bi,D)
        real(kind=double) :: x, BSCub1
        integer :: Bi, D
        BSCub1=x-(Bi-1)*h
        if ((BSCub1.LT.-2._double*h).OR.(BSCub1.GT.2._double*h)) then
            BSCub1=0._double
            return
        elseif ((BSCub1.GE.-2._double*h).AND.(BSCub1.LE.-h)) then
            BSCub1=Horner((/1._double/(6*h*h*h), 1._double/(h*h), 2._double/h, 4._double/3._double/),BSCub1,D,3)
            return
        elseif  ((BSCub1.GE.-h).AND.(BSCub1.LE.0._double)) then
            BSCub1=Horner((/-1._double/(2._double*h*h*h), -1._double/(h*h), 0._double, 2._double/3._double/),BSCub1,D,3)
            return
        elseif  ((BSCub1.GE.0._double).AND.(BSCub1.LE.h)) then
            BSCub1=Horner((/1._double/(2._double*h*h*h), -1._double/(h*h), 0._double, 2._double/3._double/),BSCub1,D,3)
            return
        elseif  ((BSCub1.GE.1._double*h).AND.(BSCub1.LE.2._double*h)) then
            BSCub1=Horner((/-1._double/(6._double*h*h*h), 1._double/(h*h), -2._double/h, 4._double/3._double/),BSCub1,D,3)
            return
        end if
    end function BSCub1

    function Horner(a,z,p,g)

        implicit none
        integer :: p, l, i, g
        real(kind=double) :: z, a(:), Horner
        real(kind=double), allocatable :: a1(:)
        allocate(a1(g+1))
        a1=a
        do l = 1, p
            do i = 2, (g + 2 - p)
                a1(i)=a1(i-1)*z+a1(i)
            end do
        end do
        Horner = a1(g+2-p)
    end function Horner


    function ProdutInterL(Bi,Bj)

        implicit none
        real(kind=double) ::  ProdutInterL, r1, r2, r3, w1, w2, w3, t1, t2, t3, IG(2), h1, alpha, betha, a, b
        integer :: ni, Cont, i, Bi, Bj
        r1 = -sqrt(3._double/5._double); r2 = 0._double;  r3 = -r1
        w1 =5._double/9._double; w2=8._double/9._double; w3=w1
        ProdutInterL=0.
        if ((m.EQ.2).and.(abs(Bi-Bj).GT.1)) return
        if ((m.EQ.4).and.(abs(Bi-Bj).GT.3)) return
        select case(m)
            case(2)
                a=max((Bi-1)*h,(Bj-1)*h,0._double)
                b=min((Bi+1)*h,(Bj+1)*h,1._double)
                ni=1
                alpha = (b-a)/2._double; betha = (b+a)/2._double
                t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                IG(1)=(w1*FProdIntSplinLin(t1,Bi,Bj)+w2*FProdIntSplinLin(t2,Bi,Bj)+w3*FProdIntSplinLin(t3,Bi,Bj))*alpha
                700 ni=2*ni; Cont = Cont+1; h1=(b-a)/ni
                IG(2)=0._double
                do i = 1, ni
                    alpha = h1/2._double; betha = ((i*h1+a)+((i-1)*h1+a))/2
                    t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                    IG(2)=IG(2)+((w1*FProdIntSplinLin(t1,Bi,Bj)+w2*FProdIntSplinLin(t2,Bi,Bj)+w3*FProdIntSplinLin(t3,Bi,Bj)))*alpha
                end do
                if (abs(IG(2)-IG(1)).GT.(10._double**(-10))) then
                    IG(1)=IG(2)
                    goto 700
                end if
                ProdutInterL = IG(2)
                return
            case(4)
                if (Bi.LE.Bj) then
                    a=max(0._double,(Bj-3)*h)
                    b=min(1._double,(Bi+1)*h)
                else
                    a=max(0._double,(Bi-3)*h)
                    b=min(1._double,(Bj+1)*h)
                end if
                ni=1
                alpha = (b-a)/2._double; betha = (b+a)/2._double
                t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                IG(1)=(w1*FProdIntLSC(t1,Bi,Bj)+w2*FProdIntLSC(t2,Bi,Bj)+w3*FProdIntLSC(t3,Bi,Bj))*alpha
                900 ni=2*ni; Cont = Cont+1; h1=(b-a)/ni
                IG(2)=0._double
                do i = 1, ni
                    alpha = h1/2._double; betha = ((i*h1+a)+((i-1)*h1+a))/2
                    t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                    IG(2)=IG(2)+((w1*FProdIntLSC(t1,Bi,Bj)+w2*FProdIntLSC(t2,Bi,Bj)+w3*FProdIntLSC(t3,Bi,Bj)))*alpha
                end do
                if (abs(IG(2)-IG(1)).GT.(10._double**(-10))) then
                    IG(1)=IG(2)
                    goto 900
                end if
                ProdutInterL = IG(2)
                return
        end select
    end function ProdutInterL


    function ProdutInterF(Bi)

       implicit none
        real(kind=double) :: ProdutInterF, r1, r2, r3, w1, w2, w3, t1, t2, t3, IG(2), h1, alpha, betha, a, b
        integer :: ni, Cont, i, Bi
        r1 = -sqrt(3._double/5._double); r2 = 0._double;  r3 = -r1
        w1 =5._double/9._double; w2=8._double/9._double; w3=w1
        select case(m)
            case(2)
                a=max((Bi-1)*h,0._double)
                b=min((Bi+1)*h,1._double)
                ni=1
                alpha = (b-a)/2._double; betha = (b+a)/2._double
                t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                IG(1)=(w1*f(t1)*BSplinLin(t1,Bi)+w2*f(t2)*BSplinLin(t2,Bi)+w3*f(t3)*BSplinLin(t3,Bi))*alpha
                7000 ni=2*ni; Cont = Cont+1; h1=(b-a)/ni
                IG(2)=0._double
                do i = 1, ni
                    alpha = h1/2._double; betha = ((i*h1+a)+((i-1)*h1+a))/2
                    t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                    IG(2)=IG(2)+(w1*f(t1)*BSplinLin(t1,Bi)+w2*f(t2)*BSplinLin(t2,Bi)+w3*f(t3)*BSplinLin(t3,Bi))*alpha
                end do
                if (abs(IG(2)-IG(1)).GT.(10._double**(-10))) then
                    IG(1)=IG(2)
                    goto 7000
                end if
                ProdutInterF = IG(2)
                return
            case(4)
                a=max(0._double,(Bi-3)*h)
                b=min(1._double,(Bi+1)*h)
                ni=1
                alpha = (b-a)/2._double; betha = (b+a)/2._double
                t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                IG(1)=(w1*f(t1)*BSCub1(t1,Bi,1)+w2*f(t2)*BSCub1(t2,Bi,1)+w3*f(t3)*BSCub1(t3,Bi,1))*alpha
                9000 ni=2*ni; Cont = Cont+1; h1=(b-a)/ni
                IG(2)=0._double
                do i = 1, ni
                    alpha = h1/2._double; betha = ((i*h1+a)+((i-1)*h1+a))/2
                    t1 = alpha*r1+betha; t2 = alpha*r2+betha; t3 = alpha*r3+betha;
                    IG(2)=IG(2)+(w1*f(t1)*BSCub1(t1,Bi,1)+w2*f(t2)*BSCub1(t2,Bi,1)+w3*f(t3)*BSCub1(t3,Bi,1))*alpha
                end do
                if (abs(IG(2)-IG(1)).GT.(10._double**(-10))) then
                    IG(1)=IG(2)
                    goto 9000
                end if
                ProdutInterF = IG(2)
                return
        end select
    end function ProdutInterF

        function ProdutInterLSC(Bi,Bj)
        real(kind=double) :: ProdutInterLSC
        real(kind=double), allocatable :: Ci(:), Cj(:)
        integer :: Bi, Bj, i, j
        integer, allocatable :: Bi1(:), Bj1(:)
        if(Bi.EQ.1) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/1, 0/)
            Ci=(/1._double, -4._double/)
        elseif(Bi.EQ.2) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/2, 0/)
            Ci=(/1._double, -1._double/)
        elseif(Bi.EQ.n+1) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/n+1, n+3/)
            Ci=(/1._double, -1._double/)
        elseif(Bi.EQ.n+2) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/n+2,n+3/)
            Ci=(/1._double, -4._double/)
        else
            allocate(Bi1(1))
            allocate(Ci(1))
            Bi1=(/Bi/)
            Ci=(/1._double/)
        end if
        if(Bj.EQ.1) then
            allocate(Bj1(2))
            allocate(Cj(2))
            Bj1=(/1, 0/)
            Cj=(/1._double, -4._double/)
        elseif(Bj.EQ.2) then
            allocate(Bj1(2))
            allocate(Cj(2))
            Bj1=(/2, 0/)
            Cj=(/1._double, -1._double/)
        elseif(Bj.EQ.n+1) then
            allocate(Bj1(2))
            allocate(Cj(2))
            Bj1=(/n+1,n+3/)
            Cj=(/1._double, -1._double/)
        elseif(Bj.EQ.n+2)then
            allocate(Bj1(2))
            allocate(Cj(2))
            Bj1=(/n+2,n+3/)
            Cj=(/1._double, -4._double/)
        else
            allocate(Bj1(1))
            allocate(Cj(1))
            Bj1=(/Bj/)
            Cj=(/1._double/)
        end if
        ProdutInterLSC=0._double
        do i = 1, size(Bi1)
            do j = 1, size(Bj1)
                ProdutInterLSC=ProdutInterLSC + Ci(i)*Cj(j) &
                                        *ProdutInterL(Bi1(i),Bj1(j))
            end do
        end do
    end function ProdutInterLSC

     function BSC(x,Bi)

        real(kind=double), allocatable :: Ci(:)
        real(kind=double) :: x, BSC
        integer :: Bi, i
        integer, allocatable :: Bi1(:)
        if(Bi.EQ.1) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/1, 0/)
            Ci=(/1._double, -4._double/)
        elseif(Bi.EQ.2) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/2, 0/)
            Ci=(/1._double, -1._double/)
        elseif(Bi.EQ.n+1) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/n+1, n+3/)
            Ci=(/1._double, -1._double/)
        elseif(Bi.EQ.n+2) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/n+2,n+3/)
            Ci=(/1._double, -4._double/)
        else
            allocate(Bi1(1))
            allocate(Ci(1))
            Bi1=(/Bi/)
            Ci=(/1._double/)
        end if
        BSC=0._double
        do i = 1, size(Bi1)
            BSC=BSC + Ci(i)*BSCub1(x,Bi1(i),1)
        end do
        deallocate(Bi1)
        deallocate(Ci)
    end function

    function ProdutInterFSC(Bi)

        real(kind=double) :: ProdutInterFSC
        real(kind=double), allocatable :: Ci(:)
        integer :: Bi, i
        integer, allocatable :: Bi1(:)
        if(Bi.EQ.1) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/1, 0/)
            Ci=(/1._double, -4._double/)
        elseif(Bi.EQ.2) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/2, 0/)
            Ci=(/1._double, -1._double/)
        elseif(Bi.EQ.n+1) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/n+1, n+3/)
            Ci=(/1._double, -1._double/)
        elseif(Bi.EQ.n+2) then
            allocate(Bi1(2))
            allocate(Ci(2))
            Bi1=(/n+2, n+3/)
            Ci=(/1._double, -4._double/)
        else
            allocate(Bi1(1))
            allocate(Ci(1))
            Bi1=(/Bi/)
            Ci=(/1._double/)
        end if
        ProdutInterFSC=0._double
        do i = 1, size(Bi1)
            ProdutInterFSC=ProdutInterFSC + Ci(i)*ProdutInterF(Bi1(i))
        end do
        deallocate(Bi1)
        deallocate(Ci)
    end function ProdutInterFSC
end program

subroutine Gauss(Maatrix,s,x,Dimen,nDiag)

    implicit none
    integer, parameter :: single = selected_real_kind(p=4)
    integer, parameter :: double = selected_real_kind(p=8)
    integer i, j, k, Diag2, nDiag, Dimen, t, t2, t3
    integer, allocatable :: t1(:)
    real(kind=double) :: Maatrix(nDiag,Dimen)
    real(kind=double) :: s(Dimen)
    real(kind=double) :: x(Dimen)
    real(kind=double) Suma, m
    Diag2 = nDiag/2
    allocate(t1(Diag2))
    do i = 2, Dimen
        do k = 1, Diag2
                t1(k)=i
        end do
        t2=Diag2-1
        t3=0
        do j =  Diag2+2, nDiag
            m=Maatrix(j,i-1)/Maatrix(Diag2+1,i-1)
            t=Diag2
            do k = j-1, j-Diag2,-1
                Maatrix(k,t1(t))=Maatrix(k,t1(t))-m*Maatrix(t,i-1)
                t=t-1
            end do
            do k = 1, t2
                t1(k)=t1(k)+1
            end do
            t2=t2-1
            s(i+t3)=s(i+t3)-s(i-1)*m
            t3=t3+1
            Maatrix(j,i-1)=0._double
        end do
    end do
    x(Dimen) = s(Dimen)/Maatrix(Diag2+1,Dimen)
    do i=Dimen-1,1,-1
        Suma=0
        do j=1, Diag2
            Suma = Suma + Maatrix(j,i)*x(i+Diag2-j+1)
        end do
        x(i)=(s(i)-Suma)/Maatrix(Diag2+1,i)
    end do
end subroutine

subroutine OUTPUTS(x,y,n,Name)
    !===================================================
    !subrutine para salvar os datos em um arquivo .dat
    !===================================================
    implicit none
    integer, parameter :: single = selected_real_kind(p=4)
    integer, parameter :: double = selected_real_kind(p=8)
    integer :: i, n
    real(kind=double) :: x(n), y(n)
    character(*) :: Name
    open(unit=1, file=Name)
    do i = 1, size(x)
        write(1,*) x(i), y(i)
    end do
    close(unit=1)
end subroutine
