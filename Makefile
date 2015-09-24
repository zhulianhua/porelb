all:  porelb.x preprocess.x postprocess.x

CC = mpiicc
INCLUDE = /opt/tecplot/include
CFLAGS = -O3 -xHost -std=c99 -g

preprocess.x: preprocess.c
	$(CC) $(CFLAGS) -o  $@  $^

porelb.x: porelb.o
	$(CC) $(CFLAGS) -o $@ $<

porelb.o: porelb.c 
	$(CC)  $(CFLAGS) -c $<

lb.o: lb.c
	$(CC)  $(CFLAGS) -c $<

postprocess.x: postprocess.c  /opt/tecplot/lib/tecio64.a
	$(CC) -I$(INCLUDE) $(CFLAGS) -o  $@  $^

touch:
	touch *.c

clean:
	rm *.o  *.x
