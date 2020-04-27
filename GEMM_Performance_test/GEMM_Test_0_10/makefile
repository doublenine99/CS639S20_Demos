CC=icpc
single_core=-Wall -mkl
gemm_test=-Wall -fopenmp -mkl=parallel -qopenmp

all: *.cpp *.h
		$(CC) -o single_core $(single_core) *.cpp
		$(CC) -o gemm_test $(gemm_test) *.cpp
single_core: *.cpp *.h
		$(CC) -o single_core $(single_core) *.cpp

gemm_test: *.cpp *.h
		$(CC) -o gemm_test $(gemm_test) *.cpp
clean:
		$(RM) gemm_test