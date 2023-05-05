CC=mpicc
CFLAGS=-lm
BUILDDIR = Build

$(BUILDDIR)/task.o:
	mkdir $(BUILDDIR)/
	$(CC) main.c func.c -o $(BUILDDIR)/task $(CFLAGS)

clean: 
	rm -r $(BUILDDIR)/