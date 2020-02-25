CC=g++
Single=-Wall
Parellel=-Wall -fopenmp


parellel: *.cpp
		$(CC) -o parellel $(Parellel) *.cpp
single: *.cpp
		$(CC) -o single $(Single) *.cpp
clean:
		$(RM) parellel single