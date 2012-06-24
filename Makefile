

COW_MPI ?= 0
CFLAGS = -Wall -g -O0 -std=c99

%.o : %.c
	$(CC) $(CFLAGS) -o $@ $< -DCOW_MPI=$(COW_MPI) -c

main : main.o cow.o
	$(CC) $(CFLAGS) -o $@ $^ -DCOW_MPI=$(COW_MPI)

clean :
	rm -rf main *.o *.dSYM
