FC=gfortran
#FC=mpif90

#FFLAGS=-fbounds-check -cpp
FFLAGS=-fbounds-check -cpp -DTIMING_P3
SDK_LIB  := $(shell xcrun --show-sdk-path 2>/dev/null)/usr/lib
LDFLAGS  = $(if $(wildcard $(SDK_LIB)),-L$(SDK_LIB),)
LIBS=

P3_SRC    = src/p3
KIN1D_SRC = src/kin1d

OBJS = microphy_p3.o routinesCld1d.o cld1d.o

clean:
	rm -f *.o *.mod fort.* execld out_p3.dat sfcprec.dat

distclean: clean
	rm -f sub

microphy_p3.o: $(P3_SRC)/microphy_p3.f90
	rm -f $@
	$(FC) -I. $(FFLAGS) -g -c -o $@ $<

routinesCld1d.o: $(KIN1D_SRC)/routinesCld1d.f90
	rm -f $@
	$(FC) -I. $(FFLAGS) -g -c -o $@ $<

cld1d.o: $(KIN1D_SRC)/cld1d.f90 routinesCld1d.o microphy_p3.o
	rm -f $@
	$(FC) -I. $(FFLAGS) -g -c -o $@ $<

execld: $(KIN1D_SRC)/execld.f90 $(OBJS)
	rm -f $@
	$(FC) -I. $(FFLAGS) -g $(LDFLAGS) -o $@ $< $(OBJS) $(LIBS)
# --> makes executable 'execld'

