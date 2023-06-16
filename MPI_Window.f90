program lb_4
	use mpi
	
	implicit none
	
	real:: t1, t2, h, tau, maxT
	integer:: steps_n, steps_tau, i, j, k, ierr, ProcNum, ProcRank, & 
			  dest, columnSize, extraColumnSize, startIndex, endIndex, mpi_status(MPI_STATUS_SIZE)
			  
	integer:: u_win, f_win
			  
	integer (KIND=MPI_ADDRESS_KIND) low_bound, extent

	real, allocatable:: u(:,:,:), u_reference(:,:,:), x(:), y(:), t(:), f(:,:,:)
	
	call MPI_INIT(ierr)
	
	call MPI_COMM_SIZE(MPI_COMM_WORLD, ProcNum, ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, ProcRank, ierr)
	
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
		write(*,*)
		write(*,*) "Введите tau: "
		read(*,*) tau
		write(*,*)
		write(*,*) "Введите h: "
		read(*,*) h
		write(*,*)
	end if
	
	call MPI_Bcast(h, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
	call MPI_Bcast(maxT, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
	call MPI_Bcast(tau, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
	
	! Количество шагов по осям x и y
	steps_n = INT(1.0/h)+1;

	! Количество шагов по оси t
	steps_tau = INT(maxT/tau)+1;
	
	! Инициализиируем индексы
	i = 1;
	j = 1;
	k = 1;
	maxT = 10;

	allocate(u(steps_n,steps_n,steps_tau), u_reference(steps_n,steps_n,steps_tau), & 
				x(steps_n), y(steps_n), t(steps_tau), f(steps_n,steps_n,steps_tau))
	

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
					MAXVAL(ABS(u_reference(:,:,:)-u(:,:,:)))
		write(*,*)
	end if
	
	call MPI_TYPE_GET_EXTENT(MPI_DOUBLE, low_bound, extent, ierr)
	
	call MPI_WIN_CREATE(u, steps_n*steps_n*steps_tau*extent, 1, MPI_INFO_NULL, MPI_COMM_WORLD, u_win, ierr)
	call MPI_WIN_CREATE(f, steps_n*steps_n*steps_tau*extent, 1, MPI_INFO_NULL, MPI_COMM_WORLD, f_win, ierr)
	
	u(:,:,:)=0.0;
	
	! Заполняем точки по осям x и y
	forall(i=1:steps_n)
		x(i)=(i-1)*h;
		y(i)=(i-1)*h;
	end forall

	! Заполняем точки по оси времени
	forall(k=1:steps_tau)
		t(k)=(k-1)*tau;
	end forall

	! Заполняем начальное граничное условие по времени
	forall(i=1:steps_n)
		u(i,j,1)=(x(i)**2+3*x(i)-9)*cos(4*y(i));
	end forall

	! Заполняем граничные условия
	forall(i=1:steps_n, k=2:steps_tau)
		u(1,i,k)=(-9)*cos(4*y(i))*exp(-2*t(k));
		u(steps_n,i,k)=(-5)*cos(4*y(i))*exp(-2*t(i));
		u(i,1,k)=(x(i)**2+3*x(i)-9)*exp(-2*t(k));
		u(i,steps_n,k)=(x(i)**2+3*x(i)-9)*cos(4.0)*exp(-2*t(k));
	end forall
	
	if (ProcRank==0) then
		t1 = MPI_WTIME();
	end if 
	
	columnSize=steps_n/ProcNum
	extraColumnSize=mod(steps_n,ProcNum)
	
	do k=2,steps_tau-1
		do i=2,steps_n-1	 
			if (ProcRank<=extraColumnSize) then
				startIndex = 1 + columnSize * (ProcRank);
				endIndex = columnSize * (ProcRank+1) + 1;
			else
				startIndex = 1 + columnSize * ProcRank;
				endIndex = steps_n;
			end if		
			
			
			if (ProcRank > 0) then				
				call MPI_WIN_FENCE(0, u_win, ierr)        
				call MPI_GET(u(i-1:i+1,startIndex:endIndex,k-1:k+1), endIndex-startIndex+1, MPI_DOUBLE, &
							 0, (endIndex-startIndex+1)*extent, endIndex-startIndex+1, &
							 MPI_DOUBLE, u_win, ierr)	
				
				
				call MPI_WIN_FREE(u_win, ierr);
				call MPI_WIN_FENCE(0, f_win, ierr)
							 
				call MPI_GET(f(i,startIndex:endIndex,k-1), endIndex-startIndex+1, MPI_DOUBLE, &
							 0, (endIndex-startIndex+1)*extent, endIndex-startIndex+1, &
							 MPI_DOUBLE, f_win, ierr)  
							
				call MPI_WIN_FREE(f_win, ierr);
			end if
			
			do j=startIndex+1,endIndex-1						
				u(i,j,k)=((u(i+1,j,k-1)-4*u(i,j,k-1)+u(i-1,j,k-1) &
							 +u(i,j+1,k-1)+u(i,j-1,k-1))/(h**2)+f(i,j,k-1)) &
							 *tau+u(i,j,k-1);	
			end do
			
			! Отправка дочерними процессами своей части матрицы
			if (ProcRank > 0) then				
				! Синхронизируем окно для u
				call MPI_WIN_FENCE(0, u_win, ierr)        
				
				call MPI_PUT(u(i-1:i+1,startIndex:endIndex,k-1:k+1), endIndex-startIndex+1, MPI_DOUBLE, &
							 0, (endIndex-startIndex+1)*extent, endIndex-startIndex+1, &
							 MPI_DOUBLE, u_win, ierr)	
				
				! Освобождаем ресурс окна
				call MPI_WIN_FREE(u_win, ierr);
			end if
		end do
	end do
	
	if (ProcRank==0) then		
		t2 = MPI_WTIME();
		
		! Выводим результат вычислений с MPI
		write(*,*) "Параллельное вычисление MPI:"

		write(*,*) "Затраченное время: ", t2-t1, " секунд"
		write(*,*) "Максимальное абсолютное отклонение от эталона", &
					MAXVAL(ABS(u_reference(:,:,:)-u(:,:,:)))
	end if 
	
	! Функция MPI_Finalize очищает все состояния, связанные с MPI.
	call MPI_Finalize(ierr)

	deallocate(u, u_reference, x, y, t, f);

end program lb_4
