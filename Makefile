


main : main.c
	gcc -Wall -g -O0 -std=c99 -o $@ $<

clean :
	rm -f main
