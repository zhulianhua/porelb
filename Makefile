all:  simgpu testpre
#APP

CC = gcc
INCLUDE = .
CFLAGS = -g -Wall
#CFLAGS = -Wall -O1

APP: testpre

testpre: testpre.o  preprocess.o lb.o
	$(CC) -I$(INCLUDE) $(CFLAGS) -o  $@  $^

testpre.o: testpre.c 
	$(CC) -I$(INCLUDE) $(CFLAGS) -c  $<

preprocess.o: preprocess.c
	$(CC) -I$(INCLUDE) $(CFLAGS) -c $<

lb.o: lb.c
	$(CC) -I$(INCLUDE) $(CFLAGS) -c $<

touch:
	touch *.c

simgpu: simgpumain.cu simgpu.cu lb.o
	#nvcc  -o  simgpu simgpumain.cu lb.o  /opt/tecplot/lib/tecio64.a -I /opt/tecplot/include  -arch=sm_13 
	nvcc  -o  simgpu simgpumain.cu lb.o  /opt/tecplot/lib/tecio64.a -I /opt/tecplot/include # -arch=sm_13 
clean:
	rm *.o testpre
