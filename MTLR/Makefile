CC=g++
#CFLAGS=-O3
CFLAGS=-g -Wall
LD=g++
#LDFLAGS=-O3
LDFLAGS=-g

all: mtlr_train mtlr_test

clean:
	rm -f mtlr_train mtlr_test *.o

mtlr_train: Sparm.o DenseVector.o SparseVector.o common.o Main.o 
	$(LD) $(LDFLAGS) Sparm.o DenseVector.o SparseVector.o common.o Main.o  -o mtlr_train

mtlr_test: Test.o Sparm.o DenseVector.o SparseVector.o common.o 
	$(LD) $(LDFLAGS) Test.o Sparm.o DenseVector.o SparseVector.o common.o  -o mtlr_test

Test.o:
	$(CC) -c $(CFLAGS) Test.cpp -o Test.o
Sparm.o:
	$(CC) -c $(CFLAGS) Sparm.cpp -o Sparm.o
DenseVector.o:
	$(CC) -c $(CFLAGS) DenseVector.cpp -o DenseVector.o
SparseVector.o:
	$(CC) -c $(CFLAGS) SparseVector.cpp -o SparseVector.o
common.o:
	$(CC) -c $(CFLAGS) common.cpp -o common.o
Main.o:
	$(CC) -c $(CFLAGS) Main.cpp -o Main.o
