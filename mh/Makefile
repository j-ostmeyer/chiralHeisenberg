FCFLAGS=-Wall -std=legacy -O3 #-flto -fPIC #-g -Og -fsanitize=address -static-libasan
FCLIBS=
FC = gfortran
MYFILES=

FCFLAGS += -Wno-unused-label -Wno-unused-dummy-argument -Wno-unused-variable

all: cheisdwf_hmc

# MAIN PROGRAM
cheisdwf_hmc: $(MYFILES:%=%.o) cheisdwf_hmc.o
	$(FC) $(FCFLAGS) -o $@ $^ $(FCLIBS)

%.o: %.f Makefile
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f cheisdwf_hmc *.o *.mod


