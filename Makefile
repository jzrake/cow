

COW_MPI ?= 0

main : main.c
	$(CC) -Wall -g -O0 -std=c99 -o $@ $< -DCOW_MPI=$(COW_MPI)

clean :
	rm -f main
