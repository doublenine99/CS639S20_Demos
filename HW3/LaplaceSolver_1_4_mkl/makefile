CC=icpc
# Single=-Wall
Laplace_1_4=-Wall -fopenmp -mkl=parallel -qopenmp


laplace_1_4: *.cpp
		$(CC) -o laplace_1_4 $(Laplace_1_4) *.cpp
clean:
		$(RM) laplace_1_4