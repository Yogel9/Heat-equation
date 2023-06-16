program lr1
    use mpi
    implicit none
    real:: h, tau, maxT, e, ac, b
    integer:: steps_n, steps_tau, i, j, k, j_k
    integer :: size, rank, ierr
    real, allocatable:: u(:,:), u_reference(:,:), x(:), t(:), f(:,:), ThreeDMatrix(:,:), matrix_DB(:,:)
    real, allocatable:: sol(:,:), vector_Db(:),  matrix_DBloc(:,:), vector_Dbloc(:),  solloc(:,:)

    double precision :: t1, t2, delta_t
    integer :: count, i_start, i_end


	
    ! Инициализация MPI
    call MPI_Init( ierr )
    call MPI_Comm_size( MPI_COMM_WORLD, size, ierr )
    call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
    
	print *, 'Процесс №', rank, 'запущен!'
	call sleep(1)
	
	call MPI_Barrier( MPI_COMM_WORLD, ierr)
	if (rank == 0) then
	print *, 'Число процессов:', size
	end if
	
	if (rank == 0) then
	write(*,*) " "
    write(*,*) "Дано:"
	write(*,*) "U(x,t)=2*x**3+5*x-9+e**5t)"
	write(*,*) "f:"
	write(*,*) "f(x,t)=12*x+5*e**5t"
	write(*,*) "НУ:"
	write(*,*) "U(x,0)=2*x**3+5*x-9"
	write(*,*) "ГУ: "
	write(*,*) "U(0,t)=-9+e**5t"
	write(*,*) "U(1,t)=-2+e**5t"
	write(*,*)
	end if
	

	if (rank == 0) then
	write(*,*) "Введите h: "
		read(*,*) h
	end if
	call MPI_Bcast(h, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    !Количество шагов по оси x 
    steps_n = INT(1.0/h)+1;
	
	
    if (rank == 0) then
    write(*,*) "Введите maxT: "
		read(*,*) maxT
	end if
	call MPI_Bcast(maxT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

	
	
	if (rank == 0) then
	write(*,*) "Введите tau: "
		read(*,*) tau
	end if
	call MPI_Bcast(tau, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	
	
	! Количество шагов по оси t
    steps_tau = INT(maxT/tau)+1;
    
    allocate (u(steps_n,steps_tau), u_reference(steps_n,steps_tau))
    allocate (x(steps_n), t(steps_tau), f(steps_n,steps_tau))
    allocate (ThreeDMatrix(steps_n,steps_n), sol(steps_n,2), solloc(steps_n,2))
    allocate (vector_Db(steps_n), matrix_DB(steps_n,steps_n), vector_Dbloc(steps_n), matrix_DBloc(steps_n,steps_n))
	
	!Заполняем шаги по сетке
	do i = 1, steps_n
		x(i) = (i - 1) * h
	end do

    do k= 1,steps_tau
        t(k)=(k-1)*tau;
    end do

	! Заполняем значения f
    do k = 1, steps_tau
		do i = 1, steps_n
        f(i,k)=-2*(EXP(-2*t(k)))*(x(i)**2-3*x(i)+10);
		end do
	end do

	! Заполняем значения эталонной функции
    do k = 1, steps_tau
		do i = 1, steps_n
        u_reference(i,k)=(EXP(-2*t(k)))*(x(i)**2-3*x(i)+9);
		end do
    end do
	
	! Заполняем начальное граничное условие по времени
    do i=1,steps_n
        u(i,1)=x(i)**2-3*x(i)+9;
    end do

	! Заполняем граничные условия
    do k=2,steps_tau
        u(1,k)=9*(EXP(-2*t(k)))
        u(steps_n,k)=7*(EXP(-2*t(k)))
    end do

	! Один поток работает все остальные ждут
	call MPI_Barrier( MPI_COMM_WORLD, ierr)
	if (rank == 0) then
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
	
	! Рассчитываем время выполнения
    delta_t = t2 - t1
	
	
	write(*,*) "Последовательное"
	write(*,*) "Затраченное время: ", delta_t, " секунд"
	write(*,*) "Отклонение от эталоной функции", &
				MAXVAL(ABS(u_reference(:,:)-u(:,:)))				
	write(*,*)
	end if
	
	call MPI_Barrier( MPI_COMM_WORLD, ierr)
	
	
	! Засекаем время перед выполнением задачи
    t1 = MPI_Wtime()
	
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
	
		!  paral
        do i=1,steps_n
            if ((i==1).OR.(i==steps_n)) then
                vector_Db(i)=(u(i,k))/ThreeDMatrix(i,i);
            else

                vector_Db(i)=(u(i,k-1)+tau*f(i,k-1))/ThreeDMatrix(i,i)
            end if
            do j=1,steps_n
                if (i/=j) then
                    matrix_DB(i,j)=-ThreeDMatrix(i,j)/ThreeDMatrix(i,i)
                else
                	! Главная диагональ - 0
                    matrix_DB(i,j)=0
                end if
            end do
        end do
		

        sol(:,2)=0.0
        sol(:,1)=vector_Db
        
        solloc(:,2)=0.0
        solloc(:,1)=vector_Db
        j_k=0
        

	call MPI_Scatter(matrix_DB, steps_n*steps_n, MPI_INTEGER, matrix_DBloc, steps_n*steps_n, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_Scatter(sol, steps_n*2, MPI_INTEGER, solloc, steps_n*2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	
	
	
	do while ((MAXVAL(ABS(solloc(:,2)-solloc(:,1)))>e).AND.(j_k<=1000))
            j_k=j_k+1;
	
            if (MOD(j_k,2)==0) then
                do i=1,steps_n
                    solloc(i,1)=dot_product(matrix_DBloc(i,:), solloc(:,2))+vector_Db(i)
                end do
            else
                do i=1,steps_n
                    solloc(i,2)=dot_product(matrix_DBloc(i,:), solloc(:,1))+vector_Db(i)
                end do
            end if
       
        end do
		 
		
	
        
       U(:,k)=sol(:,MOD(j_k,2)+1)

    end do
	
	
	call MPI_Allgather(solloc, steps_n*2, MPI_INTEGER,sol, steps_n*2, MPI_INTEGER, MPI_COMM_WORLD, ierr) 
	t2 = MPI_Wtime()	
	

    ! Рассчитываем время выполнения
    delta_t = t2 - t1
	
	
	if (rank == 0) then
	write(*,*) "Параллельно"
	write(*,*) "Затраченное время: ", delta_t, " секунд"
	write(*,*) "Отклонение от эталоной функции", &
				MAXVAL(ABS(u_reference(:,:)-u(:,:)))			
	write(*,*)
	end if
	

    

	
	! Finalize MPI
	call MPI_Finalize( ierr )
	
	write(*,*)
	deallocate(u, u_reference, x, t, f)  
	write(*,*)
end program lr1

