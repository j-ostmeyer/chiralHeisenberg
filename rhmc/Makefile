FCFLAGS=-O3
FCLIBS=
MYFILES=
MY_C_FILES=mycpu

ifneq ($(findstring login,$(HOSTNAME)),)
  #We are on DiRAC (or similar HPC facility)
  CC = icc
  FC = ifort
  #FCFLAGS += -warn all
  FCFLAGS += -heap-arrays -qopenmp -pthread #-fp-model=source
  FCFLAGS += -ipo -no-prec-div -fp-model fast=2 -xHost
else
  CC = gcc
  FC = gfortran -std=legacy
  FCFLAGS += -Wall -Wno-unused-label -Wno-unused-dummy-argument -Wno-unused-variable -Wno-surprising
  #FCFLAGS += -g -pg -fprofile-arcs -ftest-coverage #-Og #-fsanitize=address -static-libasan
  #FCLIBS += -lgcov
  #FCFLAGS += -heap-arrays # -frecursive
  FCFLAGS += -fopenmp -march=native
  #FCFLAGS += -fdefault-real-8
endif

all: cheisdwf

# MAIN PROGRAM
cheisdwf: $(MYFILES:%=%.o) $(MY_C_FILES:%=%.o) cheisdwf.o
	$(FC) $(FCFLAGS) -o $@ $^ $(FCLIBS)

%.o: %.f Makefile
	$(FC) $(FCFLAGS) -c $<
	
mycpu.o: mycpu.c Makefile
	$(CC) $(FCFLAGS) -c $<

clean:
	rm -f cheisdwf *.o *.mod


