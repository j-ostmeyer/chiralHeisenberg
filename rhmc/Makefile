FCFLAGS=-O3 #-flto -fPIC #-g -Og -fsanitize=address -static-libasan
FCLIBS=
FC = ifort
MYFILES=

ifneq ($(filter $(HOSTNAME),d3-login02 d3-login01 viz01.pri.barkla.alces.network),)
  #We are on DIaL
  #FCFLAGS += -warn all
  FCFLAGS += -heap-arrays -qopenmp -pthread #-fp-model=source
  FCFLAGS += -ipo -no-prec-div -fp-model fast=2 -xHost
else
  FC = gfortran -std=legacy
  FCFLAGS += -Wall -Wno-unused-label -Wno-unused-dummy-argument -Wno-unused-variable
  FCFLAGS += -heap-arrays # -frecursive
  FCFLAGS += -fopenmp #-march=native
  #FCFLAGS += -fdefault-real-8
endif

all: cheisdwf

# MAIN PROGRAM
cheisdwf: $(MYFILES:%=%.o) cheisdwf.o
	$(FC) $(FCFLAGS) -o $@ $^ $(FCLIBS)

%.o: %.f Makefile
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f cheisdwf *.o *.mod


