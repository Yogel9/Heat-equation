program summ
real, dimension(10000,10000):: A, B
real, dimension(10000,10000):: C
INTEGER:: i, j
real:: t1, t2

forall(i=1:10000, j=1:10000)
	A(i,j) = i+j
	B(i,j) = i+j
end forall 


call cpu_time(t1)
forall(i=1:10000, j=1:10000)
	C(i,j)=A(i,j)+B(i,j)
end forall 
call cpu_time(t2)

write(*,*) "Затраченное время: ", t2-t1, " секунд"


call cpu_time(t1)
forall(i=1:10000, j=1:10000)
	C(j,i)=A(j,i)+B(j,i)
end forall
call cpu_time(t2)

write(*,*) "Затраченное время: ", t2-t1, " секунд"
end program summ
