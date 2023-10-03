build:
	nvcc -I./inc -c cuda.cu -o cuda.o
	mpicxx -fopenmp -c main.c -o main.o
	mpicxx -fopenmp -c functions.c -o functions.o
	mpicxx -fopenmp -o finalProject  main.o functions.o cuda.o /usr/local/cuda-11.0/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./finalProject

run:
	mpiexec -np 2 ./finalProject < input.txt > output.txt

runOn2:
	mpiexec -np 2 -machinefile  mf  -map-by  node  ./finalProject