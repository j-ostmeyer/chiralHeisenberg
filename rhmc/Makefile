FCFLAGS=-Wall -std=legacy -O3 #-flto -fPIC #-g -Og -fsanitize=address -static-libasan
FCLIBS=
FC = gfortran
MYFILES=

FCFLAGS += -Wno-unused-label -Wno-unused-dummy-argument -Wno-unused-variable
FCFLAGS += -fdefault-real-8

all: cheisdwf

# MAIN PROGRAM
cheisdwf: $(MYFILES:%=%.o) cheisdwf.o
	$(FC) $(FCFLAGS) -o $@ $^ $(FCLIBS)

%.o: %.f Makefile
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f cheisdwf *.o *.mod


