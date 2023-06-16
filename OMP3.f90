program test
	use omp_lib
	implicit none
	real:: t1, t2, h, tau, maxT
	integer:: nh, ntau, i, j, k
	integer, parameter:: last_thread=7
	real, allocatable:: u(:,:,:), u_true(:,:,:), x(:), y(:), t(:), f(:,:,:)


	write(*,*) "-----------------------Initial data-----------------------"
	write(*,*) "Function U(x,y,t)=(e**(-2t))*cos(x)*cos(y)"
	write(*,*) "Function f(x,y,t)=-2*(e**(-2t))*cos(x)*cos(y)"
	write(*,*) "Initial conditions U(x,y,0)=cos(x)*cos(y)"
	write(*,*) "Boundary conditions: "
	write(*,*) "        U(0,y,t)=(e**(-2t))*cos(y)"
	write(*,*) "        U(1,y,t)=(e**(-2t))*cos(1.0)*cos(y)"
	write(*,*) "        U(x,0,t)=(e**(-2t))*cos(x)"
	write(*,*) "        U(x,1,t)=(e**(-2t))*cos(x)*cos(1.0)"
	write(*,*) "Enter (delta x and delta y) h: "
	read(*,*) h
	nh = INT(1.0/h)+1;
	write(*,*) "(Delta x and delta y) h = ", h
	write(*,*) "Enter (time maximum) maxT: "
	read(*,*) maxT
	write(*,*) "Enter (delta time) tau: "
	read(*,*) tau
	ntau = INT(maxT/tau)+1;
	allocate(u(nh,nh,ntau), u_true(nh,nh,ntau), x(nh), y(nh), t(ntau), f(nh,nh,ntau))
	

	forall(i=1:nh)
		x(i)=(i-1)*h;
		y(i)=(i-1)*h;
	end forall
	forall(k=1:ntau)
		t(k)=(k-1)*tau;
	end forall

	forall(i=1:nh, j=1:nh, k=1:ntau)
		f(i,j,k)=-2*t(k)*(EXP(-2*t(k)))*cos(x(i))*cos(y(j));
	end forall

	forall(i=1:nh, j=1:nh, k=1:ntau)
		u_true(i,j,k)=(EXP(-2*t(k)))*cos(x(i))*cos(y(j));
	end forall

	forall(i=1:nh, j=1:nh)
		u(i,j,1)=cos(x(i))*cos(y(j));
	end forall

	forall(i=1:nh, k=2:ntau)
		u(1,i,k)=(EXP(-2*t(k)))*cos(y(i));
		u(nh,i,k)=(EXP(-2*t(k)))*cos(1.0)*cos(y(i));
		u(i,1,k)=(EXP(-2*t(k)))*cos(x(i));
		u(i,nh,k)=(EXP(-2*t(k)))*cos(x(i))*cos(1.0);
	end forall

	call cpu_time(t1)
	do k=2,ntau
		do i=2,nh-1
			do j=2,nh-1
				u(i,j,k)=((u(i+1,j,k-1)-4*u(i,j,k-1)+u(i-1,j,k-1)+u(i,j+1,k-1)+u(i,j-1,k-1))/(h**2)+f(i,j,k-1))*tau+u(i,j,k-1);
			end do
		end do
	end do
	call cpu_time(t2)
		
	OPEN(10,FILE='data_3.dat')
	WRITE(10,*) h
	WRITE(10,*) nh
	WRITE(10,*) tau
	WRITE(10,*) ntau
	
	do k=1,ntau
		do i=1,nh
			do j=1,nh
				WRITE(10,*) u(i, j, k)
			end do
		end do
	end do
	CLOSE(10)

	write(*,*) "-----------------------Consistently----------------------"

	write(*,*) "Time taken by consistently calculation was ", t2-t1, " seconds"
	write(*,*) "Maximum absolute error is ", MAXVAL(ABS(u_true(:,:,:)-u(:,:,:)))
	
	u(:,:,:)=0.0;
	forall(i=1:nh, j=1:nh)
		u(i,j,1)=cos(x(i))*cos(y(j));
	end forall

	forall(i=1:nh, k=2:ntau)
		u(1,i,k)=(EXP(-2*t(k)))*cos(y(i));
		u(nh,i,k)=(EXP(-2*t(k)))*cos(1.0)*cos(y(i));
		u(i,1,k)=(EXP(-2*t(k)))*cos(x(i));
		u(i,nh,k)=(EXP(-2*t(k)))*cos(x(i))*cos(1.0);
	end forall
	
	call cpu_time(t1)
	do k=2,ntau
		!$OMP PARALLEL SHARED(u,f,tau,h)
		!$OMP WORKSHARE
		forall (i=2:nh-1,j=2:nh-1)
				u(i,j,k)=((u(i+1,j,k-1)-4*u(i,j,k-1)+u(i-1,j,k-1)+u(i,j+1,k-1)+u(i,j-1,k-1))/(h**2)+f(i,j,k-1))*tau+u(i,j,k-1);
			
		end forall
		!$OMP END WORKSHARE
		!$OMP END PARALLEL
	end do
	call cpu_time(t2)

	write(*,*) "-----------------------Paralelly----------------------"

	write(*,*) "Time taken by consistently calculation was ", t2-t1, " seconds"
	write(*,*) "Maximum absolute error is ", MAXVAL(ABS(u_true(:,:,:)-u(:,:,:)))

	deallocate(u, u_true, x, y, t, f)

end program test

