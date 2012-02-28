# PREDEFINED VARIABLES
CFLAGS = -c -I. -O3 -Wall -Wno-unused -Wno-unused-result
LFLAGS = -lm
CC = mpicc

# GLOBAL VARIABLES
INITFILE = init.txt
INPUTFILE = entrada.dat
OUTPUTFILE = salida

# NUMBER OF PROCESS
NPROC = 2

# MODULES
MODULES = lib_simulacion integradores
MODULESC = $(MODULES:=.c)
MODULESH = $(MODULES:=.h)
MODULESO = $(MODULES:=.o)
PROGRAM = simulacion

# LINKING
$(PROGRAM).out: $(PROGRAM).o $(MODULESO)
	$(CC) $^ -o $@ $(LFLAGS)
	rm -rf $(MODULESO) $(PROGRAM).o

# MODULES' RULES
$(MODULES): $(MODULESC) $(MODULESH)
	$(CC) $(CFLAGS) $@.c -o $@.o

# PROGRAM COMPILATION
$(PROGRAM).o: $(PROGRAM).c $(MODULESH)
	$(CC) $(CFLAGS) $< -o $@

run:
	mpiexec -np $(NPROC) ./$(PROGRAM).out $(INITFILE) $(INPUTFILE) $(OUTPUTFILE) 

clean:
	rm -rf *~ *.out *.o *#

all: $(PROGRAM).out run


