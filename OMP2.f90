!  MatrixMultipl.f90 
!
!  FUNCTIONS:
!  MatrixMultipl - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MatrixMultipl
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program MatrixMultipl
    
    USE OMP_LIB

    implicit none

    ! Variables

    integer :: n, i, j, k, loop, ID, phase

    real, allocatable :: A(:,:), Bseq(:,:), Bpar(:,:)
    real :: left, right, top, bottom
	logical :: test
    
	integer(kind = OMP_lock_kind), dimension(5) :: lcks

    real :: start_time, finish_time
    
    write(*,*) "Enter n"
    read(*,*) n
    
    allocate (A(n,n))
    allocate (Bseq(n,n))
    allocate (Bpar(n,n))
	
	call OMP_SET_NUM_THREADS(4)
	call OMP_SET_NESTED(.TRUE.)
	
	write(*,*) "Max Threads = ", OMP_GET_MAX_THREADS()
    
    call random_seed

    CALL RANDOM_NUMBER(A)
	
	!****************************************************************************
    !
    ! Последовательно
    !
    !****************************************************************************
    
	Bseq = 0.0
		
    call cpu_time(start_time)
	
	! Сглаживаем несколько раз
	do phase= 1,5
		do i=1,n 
			do j=1,n
				if (i==1) then 
					top = A(n,j)
				else 
					top = A(i-1, j)
				end if
				if (i==n) then 
					bottom = A(1,j)
				else 
					bottom = A(i+1,j)
				end if
				if (j==1) then 
					left = A(i,n)
				else
					left = A(i,j-1)
				end if
				if (j==n) then 
					right = A(i,n)
				else 
					right = A(i, j+1)    
				end if
				Bseq(i,j)=(A(i,j)+left+right+top+bottom)/5
			end do 
		end do
	end do
    
    call cpu_time(finish_time)
	
	print *, 'Sequential calculation'

    print *, 'Calculation time = ', finish_time - start_time
	
	print *, ''
	
	!****************************************************************************
    !
    ! Параллельно
    !
    !****************************************************************************
	
	do i=1,5 
		call OMP_INIT_LOCK(lcks(i))
	end do
	
	Bpar = 0.0
		
    call cpu_time(start_time)
	
	! Сглаживаем несколько раз
	do phase= 1,5
		!$OMP PARALLEL SHARED(lcks) PRIVATE(top, bottom, left, right, ID, test)
		top = 0
		bottom = 0
		left = 0
		right = 0
		!$OMP DO
		do loop=1, OMP_GET_MAX_THREADS()
			ID = OMP_GET_THREAD_NUM() + 1
				
			CALL OMP_SET_LOCK(lcks(ID))
			!print *, 'My thread id is', ID
			CALL OMP_UNSET_LOCK(lcks(ID))
			
			do while (.NOT. OMP_TEST_LOCK(lcks(ID)))
				! We do not yet have the lock
				! so we must do something else
				print *, 'Waiting...', ID
			end do
			
			do i=1+((ID-1)*n/4),ID*n/4				
				do j=1,n
					if (i==1) then 
						top = A(n,j)
					else 
						top = A(i-1, j)
					end if
					if (i==n) then 
						bottom = A(1,j)
					else 
						bottom = A(i+1,j)
					end if
					if (j==1) then 
						left = A(i,n)
					else
						left = A(i,j-1)
					end if
					if (j==n) then 
						right = A(i,n)
					else 
						right = A(i, j+1)    
					end if
					
					Bpar(i,j)=(A(i,j)+left+right+top+bottom)/5	
				end do
			end do 
			
			CALL OMP_UNSET_LOCK(lcks(ID))
			
		end do
		!$OMP END DO
		!$OMP END PARALLEL
	end do
    
    call cpu_time(finish_time)
	
	do i=1,5 
		call OMP_DESTROY_LOCK(lcks(i));
	end do
	
	print *, 'Parallel calculation'

    print *, 'Calculation time = ', finish_time - start_time
	
	
	print *, ''
	
	print *, 'maxval(abs(Bseq-Bpar)) = ', maxval(abs(Bseq-Bpar))

    write(*,*) "Press num button..."
    read(*,*) n
    

    end program MatrixMultipl
