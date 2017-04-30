OBJS = main.o ode_method.o
CC= g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
main : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o main
main.o: main.cpp ode_method.h
	$(CC) $(CFLAGS) main.cpp
ode_method.o: ode_method.cpp ode_method.h
	$(CC) $(CFLAGS) ode_method.cpp

clean:
	\rm *.o main
