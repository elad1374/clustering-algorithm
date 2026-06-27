
all: symnmf.o symnmf.h prototypes.h 
	gcc -o symnmf symnmf.o -lm

symnmf.o: symnmf.c
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c symnmf.c -lm
