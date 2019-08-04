bin = $(HOME)/bin
CC_PLUS = g++  
CC = gcc  
C_FLAGS = -O2

all: $(bin)/nucleation_tracker

$(bin)/nucleation_tracker:  cmdline.c nucleation_tracker.cpp
	$(CC) $(C_FLAGS) -c cmdline.c
	$(CC_PLUS) $(C_FLAGS) -c nucleation_tracker.cpp
	$(CC_PLUS) -o $(bin)/nucleation_tracker cmdline.o nucleation_tracker.o

clean: 
	rm *~
	rm *.o
