program lb_4
    use omp_lib
    implicit none
    real:: t1, t2, h, tau, maxT, e, ac, b
    integer:: steps_n, steps_tau, i, j, k, j_k
    
    real, allocatable:: u(:,:), u_reference(:,:), x(:), t(:), f(:,:), & 
    ThreeDMatrix(:,:), matrix_DB(:,:), sol(:,:), vector_Db(:)

    
    write(*,*) "Дано:"
	write(*,*) "U(x,t)=e**(-2t)*(x**2-3*x+9)"
	write(*,*)
	
	write(*,*) "Находим f:"
	write(*,*) "f(x,y,t)=-2*e**(-2t)*(x**2-3*x+10)"
	write(*,*)
	
	write(*,*) "Начальные условия:"
	write(*,*) "U(x,0)=x**2-3*x+9"
	write(*,*)
	
	write(*,*) "Граничные условия: "
	write(*,*) "U(0,t)=9*e**(-2t)"
	write(*,*) "U(1,t)=7*e**(-2t)"
	write(*,*)
	
	write(*,*) "Введите h: "
	read(*,*) h
	write(*,*)
    
    !Количество шагов по оси x 
    steps_n = INT(1.0/h)+1;
    
    write(*,*) "Введите ограничение по времени maxT: "
	read(*,*) maxT
	write(*,*)

	write(*,*) "Введите tau: "
	read(*,*) tau
	write(*,*)
	
	! Количество шагов по оси t
    steps_tau = INT(maxT/tau)+1;
    
    allocate(u(steps_n,steps_tau), u_reference(steps_n,steps_tau), & 
    x(steps_n), t(steps_tau), f(steps_n,steps_tau), & 
    ThreeDMatrix(steps_n,steps_n), & 
    sol(steps_n,2), vector_Db(steps_n), & 
    matrix_DB(steps_n,steps_n))
    
	! Заполняем точки по оси x
    forall(i=1:steps_n)
        x(i)=(i-1)*h;
    end forall
    
    ! Заполняем точки по оси времени
    forall(k=1:steps_tau)
        t(k)=(k-1)*tau;
    end forall

	! Заполняем значения f
    forall(i=1:steps_n, k=1:steps_tau)
        f(i,k)=-2*(EXP(-2*t(k)))*(x(i)**2-3*x(i)+10);
    end forall

	! Заполняем значения эталонной функции
    forall(i=1:steps_n, k=1:steps_tau)
        u_reference(i,k)=(EXP(-2*t(k)))*(x(i)**2-3*x(i)+9);
    end forall
    
	! Заполняем начальное граничное условие по времени
    forall(i=1:steps_n)
        u(i,1)=x(i)**2-3*x(i)+9;
    end forall

	! Заполняем граничные условия
    forall(k=2:steps_tau)
        u(1,k)=9*(EXP(-2*t(k)))
        u(steps_n,k)=7*(EXP(-2*t(k)))
    end forall

    call cpu_time(t1)

    do k=2,steps_tau
        
        ThreeDMatrix=0.0;
        ThreeDMatrix(1,1)=1.0;
        ThreeDMatrix(steps_n, steps_n)=1.0;

		! Главная диагональ заполняется 1+2*tau/(h**2), соседние - -tau/(h**2) (Трёхдиагональная матрица с U)
        do i=2,steps_n-1
            ThreeDMatrix(i,i-1)=-tau/(h**2)
            ThreeDMatrix(i,i)=1+2*tau/(h**2)
            ThreeDMatrix(i,i+1)=-tau/(h**2);
        end do

        do i=1,steps_n
        	! Если первый или последний элемент, то u(i,k)/главную диагональ ThreeDMatrix
        	! Аналог вектора D*b (Метод Якоби)
            if ((i==1).OR.(i==steps_n)) then
                vector_Db(i)=(u(i,k))/ThreeDMatrix(i,i);
            else
            	! Иначе - u(i,k-1) + tau*f(i,k-1)) /главную диагональ ThreeDMatrix
                vector_Db(i)=(u(i,k-1)+tau*f(i,k-1))/ThreeDMatrix(i,i)
            end if
            do j=1,steps_n
                if (i/=j) then
                	! Еслси элемент не на главной диагонали, то -ThreeDMatrix(i,j)/главную диагональ ThreeDMatrix
                	! Аналог матрицы D*B (метод Якоби)
                    matrix_DB(i,j)=-ThreeDMatrix(i,j)/ThreeDMatrix(i,i)
                else
                	! Главная диагональ - 0
                    matrix_DB(i,j)=0
                end if
            end do
        end do
		
		! Для свопа (два соседних решения)
        sol(:,2)=0.0
        sol(:,1)=vector_Db
        j_k=0
        
        do while ((MAXVAL(ABS(sol(:,2)-sol(:,1)))>e).AND.(j_k<=1000))
            j_k=j_k+1;
            if (MOD(j_k,2)==0) then
                do i=1,steps_n
                    sol(i,1)=dot_product(matrix_DB(i,:), sol(:,2))+vector_Db(i)
                end do
            else
                do i=1,steps_n
                    sol(i,2)=dot_product(matrix_DB(i,:), sol(:,1))+vector_Db(i)
                end do
            end if
        end do

        U(:,k)=sol(:,MOD(j_k,2)+1)

    end do
    call cpu_time(t2)
        
    OPEN(10,FILE='lb_4_plot_data.dat')
    WRITE(10,*) h
    WRITE(10,*) steps_n
    WRITE(10,*) tau
    WRITE(10,*) steps_tau
    
    do k=1,steps_tau
        do i=1,steps_n
           WRITE(10,*) u(i, k)
        end do
    end do
    CLOSE(10)

	write(*,*) "Последовательное вычисление:"

	write(*,*) "Затраченное время: ", t2-t1, " секунд"
	write(*,*) "Максимальное абсолютное отклонение от эталона", &
				MAXVAL(ABS(u_reference(:,:)-u(:,:)))
	write(*,*)

    u=0.0
    
    forall(i=1:steps_n)
        u(i,1)=x(i)**2-3*x(i)+9;
    end forall

    forall(k=2:steps_tau)
        u(1,k)=9*(EXP(-2*t(k)))
        u(steps_n,k)=7*(EXP(-2*t(k)))
    end forall

    call cpu_time(t1)

    do k=2,steps_tau
        
        ThreeDMatrix=0.0;ThreeDMatrix(1,1)=1.0;ThreeDMatrix(steps_n, steps_n)=1.0;

        do i=2,steps_n-1
            ThreeDMatrix(i,i-1)=-tau/(h**2)
            ThreeDMatrix(i,i)=1+2*tau/(h**2);
            ThreeDMatrix(i,i+1)=-tau/(h**2)
        end do

        do i=1,steps_n
            if ((i==1).OR.(i==steps_n)) then
                vector_Db(i)=(u(i,k))/ThreeDMatrix(i,i);
            else
                vector_Db(i)=(u(i,k-1)+tau*f(i,k-1))/ThreeDMatrix(i,i);
            end if
            do j=1,steps_n
                if (i/=j) then
                    matrix_DB(i,j)=-ThreeDMatrix(i,j)/ThreeDMatrix(i,i);
                else
                    matrix_DB(i,j)=0;
                end if
            end do
        end do

        sol(:,2)=0.0;sol(:,1)=vector_Db;j_k=0;
        do while ((MAXVAL(ABS(sol(:,2)-sol(:,1)))>e).AND.(j_k<=1000))
            j_k=j_k+1;
            if (MOD(j_k,2)==0) then
                !$OMP PARALLEL DO SHARED(vector_Db, sol, matrix_DB)
                do i=1,steps_n
                    sol(i,1)=dot_product(matrix_DB(i,:), sol(:,2))+vector_Db(i)
                end do
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO SHARED(vector_Db, sol, matrix_DB)
                do i=1,steps_n
                    sol(i,2)=dot_product(matrix_DB(i,:), sol(:,1))+vector_Db(i)
                end do
                !$OMP END PARALLEL DO
            end if
        end do

        U(:,k)=sol(:,MOD(j_k,2)+1)

    end do
    call cpu_time(t2)

	write(*,*) "Параллельное вычисление:"

	write(*,*) "Затраченное время: ", t2-t1, " секунд"
	write(*,*) "Максимальное абсолютное отклонение от эталона", &
				MAXVAL(ABS(u_reference(:,:)-u(:,:)))

    deallocate(u, u_reference, x, t, f, ThreeDMatrix, sol, vector_Db, matrix_DB)

end program lb_4

