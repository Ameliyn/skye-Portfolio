CC=cc
EXEC=bitfun
OBJS=driver.o bitlib.o
OPTS=-std=c99 -Wall
SOURCES=*.c

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) -o $(EXEC)

.c.o:
	$(CC) $(OPTS) -c *.c

clean:
	rm -f *.o $(EXEC) a.out
