program LR2
	use mpi
	
	
	real:: t1, t2, h, tau, maxT
	integer:: steps_n, steps_tau, i, j, k, IErr, ProcNum, ProcRank, & 
			  dest, columnSize, extraColumnSize, startIndex, endIndex, MPI_STATUS(MPI_STATUS_SIZE)

	real, allocatable:: u(:,:,:), u_reference(:,:,:), x(:), y(:), t(:), f(:,:,:), vector(:), res(:), calc_matrix(:,:,:)
	
	! Тип сообщения родительскому процессу
	integer, parameter::MASTER=1
	! Тип сообщения между рабочими процессами
	integer, parameter::WORKER_F=2	
	! Тип сообщения между рабочими процессами
	integer, parameter::WORKER_L=3
	
	call MPI_INIT(IErr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, ProcNum, IErr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, ProcRank, IErr)
	
	! Выводим информацию
	! Вводим данные
	if (ProcRank==0) then
		write(*,*) "Процессов в работе:", ProcNum
		write(*,*) "Дано:"
		write(*,*) "U(x,y,t)=5*x**2+sin(y)*exp(3*t)"
		write(*,*) "f:"
		write(*,*) "f(x,y,t)=3*sin(y)*exp(3*t)+exp(3*t)*sin(y)-10"
		write(*,*) "Начальные условия:"
		write(*,*) "U(x,y,0)=5*x**2+sin(y)"
		write(*,*) "Граничные условия: "
		write(*,*) "U(0,y,t)=sin(y)*exp(3*t)"
		write(*,*) "U(1,y,t)=5+sin(y)*exp(3*t)"
		write(*,*) "U(x,0,t)=5*x**2"
		write(*,*) "U(x,1,t)=5*x**2+sin(1.0)*exp(3*t)"
		write(*,*) "Введите tau: "
		read(*,*) tau
		write(*,*)
		write(*,*) "Введите h: "
		read(*,*) h
		write(*,*)
	end if
	
	i = 1;
	j = 1;
	k = 1;
	maxT = 10;

	call MPI_Bcast(h, 1, MPI_INT, 0, MPI_COMM_WORLD, IErr)
	call MPI_Bcast(maxT, 1, MPI_INT, 0, MPI_COMM_WORLD, IErr)
	call MPI_Bcast(tau, 1, MPI_INT, 0, MPI_COMM_WORLD, IErr)
	
	! Количество шагов по осям x и y
	steps_n = INT(1.0/h)+1;

	! Количество шагов по оси t
	steps_tau = INT(maxT/tau)+1;
	
	allocate(u(steps_n,steps_n,steps_tau), u_reference(steps_n,steps_n,steps_tau), & 
				x(0:steps_n+1), y(steps_n), t(steps_tau), f(steps_n,steps_n,steps_tau))
	
	if (ProcRank==0) then
		! Заполняем точки по осям x и y
		forall(i=1:steps_n)
			x(i)=(i-1)*h;
			y(i)=(i-1)*h;
		end forall

		! Заполняем точки по оси времени
		forall(k=1:steps_tau)
			t(k)=(k-1)*tau;
		end forall

		! Заполняем значения f и u_ref
		forall(j=1:steps_n, i=1:steps_n, k=1:steps_tau)
			f(j,i,k)=3*sin(y(i))*exp(3*t(k))+exp(3*t(k))*sin(y(i))-10;
			u_reference(i,j,k)=5*x(i)**2+sin(y(i))*exp(3*t(i));
		end forall

		! Заполняем начальное граничное условие по времени
		forall(i=1:steps_n)
			u(i,j,1)=(x(i)**2+3*x(i)-9)*cos(4*y(i));
		end forall

		! Заполняем граничные условия
		forall(i=1:steps_n, k=2:steps_tau)
			u(1,i,k)=sin(y(i))*exp(3*t(k));
			u(steps_n,i,k)=5+sin(y(i))*exp(3*t(k));
			u(i,1,k)=5*x(i)**2;
			u(i,steps_n,k)=5*x(i)**2+sin(1.0)*exp(3*t(k));
		end forall
	
		call cpu_time(t1)
		
		do k=2,steps_tau
			do i=2,steps_n-1
				do j=2,steps_n-1
					! Заполняем оставшиеся значения
					u(i,j,k)=((u(i+1,j,k-1)-4*u(i,j,k-1)+u(i-1,j,k-1)+u(i,j+1,k-1)+u(i,j-1,k-1))/(h**2)+f(i,j,k-1))*tau+u(i,j,k-1);
				end do
			end do
		end do
		
		call cpu_time(t2)

		write(*,*) "Последовательное вычисление:"

		write(*,*) "Затраченное время: ", t2-t1, " секунд"
		write(*,*) "Максимальное абсолютное отклонение от эталона", &
					MAXVAL(ABS(u_reference(steps_n,:,:)-u(steps_n,:,:)))
		write(*,*)
	end if
	
	u(:,:,:)=0.0;
	
	columnSize=steps_n/ProcNum
	extraColumnSize=mod(steps_n,ProcNum)
	
	allocate(vector(columnSize), res(0:columnSize+1), calc_matrix(3,0:columnSize+1,3))
	

	t1 = MPI_WTIME();

	do k=2,steps_tau-1
		do i=2,steps_n-1	
			if (ProcRank==0) then
				do dest=1, ProcNum-1
					! Вычисляем индексы
					if (dest<=extraColumnSize) then
						startIndex = 1 + columnSize * (dest - 1);
						endIndex = columnSize * dest + 1;
					else
						startIndex = 1 + columnSize * (dest - 1);
						endIndex = steps_n;
					end if
					
					call MPI_

					Recv(res, columnSize + 2, MPI_DOUBLE, dest, MASTER, MPI_COMM_WORLD, MPI_STATUS, IErr)
						
					u(1,startIndex:endIndex,k)=res(1:columnSize);
				end do
			end if 
			
			if (ProcRank>0) then	
				! Вычисляем индексы
				if (ProcRank<=extraColumnSize) then
					startIndex = 1 + columnSize * (ProcRank - 1);
					endIndex = columnSize * ProcRank + 1;
				else
					startIndex = 1 + columnSize * (ProcRank - 1);
					endIndex = steps_n;
				end if			
				
				! Инициализиируем
				calc_matrix(:,:,:)=0.0
				
				if (startIndex > 1 .AND. endIndex < steps_n) then
					calc_matrix(:,0:columnSize+1,:) = u(1:i+1,startIndex-1:endIndex+1,k-1:k+1);
				else if (startIndex > 1 .AND. endIndex == steps_n) then
					calc_matrix(:,0:columnSize,:) = u(1:i+1,startIndex-1:endIndex,k-1:k+1);
				else if (startIndex == 1 .AND. endIndex < steps_n) then
					calc_matrix(:,1:columnSize+1,:) = u(1:i+1,startIndex:endIndex+1,k-1:k+1);
				else
					calc_matrix(:,1:columnSize,:) = u(1:i+1,startIndex:endIndex,k-1:k+1);
				end if
				
				vector(:) = f(k,startIndex:endIndex,k-1);
				
				res(:) = 0.0;
				res(1:columnSize) = calc_matrix(2,startIndex:endIndex,2);
				
				do j=2,columnSize-1						
					res(j)=((calc_matrix(3,j,1)-4*calc_matrix(2,j,1) &
							+calc_matrix(1,j,1)+calc_matrix(2,j+1,1) &
							+calc_matrix(2,j-1,1))/(h**2)+vector(j)) &
							*tau+calc_matrix(2,j,1);	
				end do
				
				!  Отправляем Столбец 1 к ID-1
				if ( ProcRank > 1 ) then
					call MPI_Send ( res(1), 1, MPI_DOUBLE, ProcRank-1, WORKER_F, MPI_COMM_WORLD, IErr)
				end if
				
				!  Получаем Столбец N+1 от ID+1
				if ( ProcRank < ProcNum - 1 ) then
					call MPI_Recv ( res(columnSize+1), 1, MPI_DOUBLE, ProcRank+1, WORKER_F, MPI_COMM_WORLD, MPI_STATUS, IErr)
					
					j=columnSize;
					res(j)=((calc_matrix(3,j,1)-4*calc_matrix(2,j,1) &
							+calc_matrix(1,j,1)+calc_matrix(2,j+1,1) &
							+calc_matrix(2,j-1,1))/(h**2)+vector(j)) &
							*tau+calc_matrix(2,j,1);
				else
					j=columnSize;
					res(j)=((calc_matrix(3,j,1)-4*calc_matrix(2,j,1) &
							+calc_matrix(1,j,1)+calc_matrix(2,j+1,1) &
							+calc_matrix(2,j-1,1))/(h**2)+vector(j)) &
							*tau+calc_matrix(2,j,1);
				end if
				
				!  Отправляем Столбец N к ID+1
				if ( ProcRank < ProcNum - 1 ) then
					call MPI_Send ( res(columnSize), 1, MPI_DOUBLE, ProcRank+1, WORKER_L, MPI_COMM_WORLD, IErr)
				end if
	
				!  Получаем Столбец 0 from ID-1.
				if ( ProcRank > 1 ) then
					call MPI_Recv (res(0), 1,  MPI_DOUBLE, ProcRank-1, WORKER_L, MPI_COMM_WORLD, MPI_STATUS, IErr)
					
					j=1;
					res(j)=((calc_matrix(3,j,1)-4*calc_matrix(2,j,1) &
							+calc_matrix(1,j,1)+calc_matrix(2,j+1,1) &
							+calc_matrix(2,j-1,1))/(h**2)+vector(j)) &
							*tau+calc_matrix(2,j,1);
				else
					j=1;
					res(j)=((calc_matrix(3,j,1)-4*calc_matrix(2,j,1) &
							+calc_matrix(1,j,1)+calc_matrix(2,j+1,1) &
							+calc_matrix(2,j-1,1))/(h**2)+vector(j)) &
							*tau+calc_matrix(2,j,1);
				end if
				
				! Отправляем результат u
				call MPI_Send(res(:), columnSize+2, MPI_DOUBLE, 0, MASTER, MPI_COMM_WORLD, IErr)			
			end if
		end do
	end do
	t2 = MPI_WTIME();

	if (ProcRank==0) then		
		write(*,*) "Параллельное вычисление MPI:"

		write(*,*) "Затраченное время: ", t2-t1, " секунд"
		write(*,*) "Максимальное абсолютное отклонение от эталона", &
					MAXVAL(ABS(u_reference(steps_n,:,:)-u(steps_n,:,:)))
	end if 
	
	! Функция MPI_Finalize очищает все состояния, связанные с MPI.
	call MPI_Finalize(IErr)

	deallocate(u, u_reference, x, y, t, f, vector, res, calc_matrix);

end program LR2
