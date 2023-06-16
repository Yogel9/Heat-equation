
    program MatrixMultipl
    
    USE OMP_LIB

    implicit none

    ! Variables

    integer :: n, i, j, k

    real, allocatable :: A(:,:), D(:,:), difference_DA(:,:), b(:), x0(:), xsol(:), &
    vector_difference(:), matrix_B(:,:), g(:), iterative_x(:,:)
    
    real :: summ, norma, e, l

    real :: start_time, finish_time
    
    write(*,*) "Enter n"
    read(*,*) n
    
    allocate (A(n,n))
    allocate (D(n,n))
	allocate (difference_DA(n,n))
    allocate (iterative_x(2,n))
    allocate (b(n))
    allocate (x0(n))
    allocate (xsol(n))
    allocate (vector_difference(n))
    allocate (matrix_B(n,n))
    allocate (g(n))
    
    
    print *, 'Start program'
    
   !формируем матрицу A
    CALL random_seed
    CALL RANDOM_NUMBER(A)
    
    
    A=A-0.5
    
    do i=1,n
    	summ=0
    	do j=1,n
        	if (i.ne.j) then
        	    summ=summ+abs(A(i,j))
        	end if
    	end do
	    A(i,i)=summ
    end do
    
    
    !формируем решение xsol
    CALL RANDOM_NUMBER(xsol)
	!находим b
    call multiplication_matrix_vector(A,xsol,b)
    
	!приводим преобразование системы уравнений к итерационному виду
    D=0.0
    
    do i=1,n
	    D(i,i)=A(i,i)!матрица D (главная диагональ матрицы А)
    end do

    difference_DA=D-A
    
    do i=1,n
        g(i) = b(i) / D(i,i)
    end do
    
    do i=1,n
        do j=1,n 
            matrix_B(i,j) = difference_DA(i,j) / D(i,i)
        end do
    end do
    
    CALL RANDOM_NUMBER(x0)
	
	
    ! Последовательно
   
    i = 1
    j = 2 
    k=0
    
    !Порядок нормы
    l=2
    
    !Точность
    e=0.01
    
    !Начальная норма
    norma=1
    
    iterative_x(i,1:n)=x0
    iterative_x(j,1:n)=0.0
    
    call cpu_time(start_time)
    
    do while (norma>=e)
    
        call multiplication_matrix_vector(matrix_B,iterative_x(i,1:n),iterative_x(j,1:n))
        iterative_x(j,1:n) = iterative_x(j,1:n)+g
        vector_difference=iterative_x(j,1:n)-iterative_x(i,1:n)
        call find_norma(vector_difference, n, l, norma)
        
        call i_swap(i,j)
        
    end do
    
    call cpu_time(finish_time)
    

    print *, 'Sequential calculation'
    print *, 'Calculation time = ', finish_time - start_time
    print *, ''
    
   
    ! Параллельно
    
    call cpu_time(start_time)
    
    i = 1
    j = 2 
    k=0
    
    !Порядок нормы
    l=2
    
    !Точность
    e=0.01
    
    !Начальная норма
    norma=1
    
    iterative_x(i,1:n)=x0
    iterative_x(j,1:n)=0.0
    
    do while (norma>=e)
    
        call multiplication_matrix_vector_omp(matrix_B,iterative_x(i,1:n),iterative_x(j,1:n))
        iterative_x(j,1:n) = iterative_x(j,1:n)+g
        vector_difference=iterative_x(j,1:n)-iterative_x(i,1:n)
        call find_norma_omp(vector_difference, n, l, norma)
        
        call i_swap(i,j)
        
    end do
    
    call cpu_time(finish_time)
    
    print *, 'Parallel calculation'
    print *, 'Calculation time = ', finish_time - start_time
    

    write(*,*) "Write something numb..."
    read(*,*) n
    
contains
    
    subroutine multiplication_matrix_vector(matrix, vector, result_vector)
        integer :: i, length
        real, dimension(:,:), intent(in) :: matrix
        real, dimension(:), intent(in) :: vector
        real, dimension(:), intent(out) :: result_vector
        
		length=size(matrix,1)
        result_vector=0.0
        
    	do i=1,length
    		result_vector(i) = dot_product(matrix(i,1:length),vector)
    	end do

    end subroutine multiplication_matrix_vector 
    
    subroutine multiplication_matrix_vector_omp(matrix, vector, result_vector)
        use OMP_LIB

        integer :: i, j, length_m, length_v
        real summ
        real, dimension(:,:), intent(in) :: matrix
        real, dimension(:), intent(in) :: vector
        real, dimension(:), intent(out) :: result_vector
        
        length_m=size(matrix,1)
        length_v=size(vector,1)
        
        result_vector=0.0
        
        !$OMP PARALLEL DO SHARED(result_vector, matrix)
    	do i=1,length_m
    		result_vector(i) = dot_product(matrix(i,1:length_m),vector)
    	end do
    	!$OMP END PARALLEL DO

    end subroutine multiplication_matrix_vector_omp  
    
    subroutine find_norma(vector, n, l, norma)
        integer, intent(in):: n
        integer :: i, j
        real, intent(in) :: l
        real, dimension(:), intent(in) :: vector
        real :: summ
        real, intent(out) :: norma
        
    	summ=0;
                
        do i=1,n
            summ=summ+vector(i)**l
        end do
        
        norma=summ**(1/l) 
        
        return
    end subroutine find_norma
        
    subroutine find_norma_omp(vector, n, l, norma)
        use OMP_LIB
    
        integer, intent(in):: n
        integer :: i, j
        real, intent(in) :: l
        real, dimension(:), intent(in) :: vector
        real :: summ
        real, intent(out) :: norma
        
    	summ=0;
                
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(summ, i)
        do i=1,n
            summ=summ+vector(i)**l
        end do
        !$OMP END PARALLEL DO
        
        norma=summ**(1/l) 
        
        return
    end subroutine find_norma_omp
    
    subroutine i_swap(lhs,rhs)
      integer, intent(inout) :: lhs,rhs
      integer                :: temp
      temp = lhs; lhs = rhs; rhs = temp
   end subroutine i_swap

    end program MatrixMultipl
