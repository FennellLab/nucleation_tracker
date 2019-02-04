bin = $(HOME)/bin
CC_PLUS = g++  
C_FLAGS = -O2

all: $(bin)/nucleation_tracker

$(bin)/nucleation_tracker:  nucleation_tracker.cpp
	$(CC_PLUS) nucleation_tracker.cpp $(C_FLAGS) -o $(bin)/nucleation_tracker

clean: 
	rm *~
	rm *.o
