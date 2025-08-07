FC=gfortran
#FC=mpif90
#FC=f90
#FC=pgf90

#FFLAGS=-C -g
##FFLAGS=-C -g -qflttrap=invalid:nanq:ov:zerodivide:enable:imp:overflow  #IBM
##FFLAGS=-fdefault-real-8 -Wno-line-truncation -fbounds-check -finit-real=nan 
##FFLAGS=-fdefault-real-8 -Wno-line-truncation 
##FFLAGS=-fdefault-real-8 -Wno-line-truncation -fbounds-check -finit-real=nan -Wall -Wno-conversion -Wno-unused-dummy-argument -std=f2008 -Wunused-parameter
#FFLAGS=-fdefault-real-8 -Wno-line-truncation -fbounds-check -finit-real=nan -Wall -Wno-unused-dummy-argument -std=f2008 -Wno-intrinsics-std -Wunused-parameter

#FFLAGS=-fbounds-check -cpp 
FFLAGS=-fbounds-check -cpp -DTIMING_P3
LDFLAGS=
LIBS=



%.o : %.mod
.SUFFIXES: .f90 .o

.f90.o:
	rm -f $@
	$(FC) -I. -I. $(FFLAGS) -g -c -o $@ $<
.o:
	rm -f $@
	$(FC) -I. -I. $(FFLAGS) -g  $(LDFLAGS) -o $@ $< $(OBJS) $(LIBS)
default:

clean:
	rm -f *.o *.mod fort.* execld out_p3.dat sfcprec.dat

OBJS = microphy_p3.o routinesCld1d.o   cld1d.o
exec:
	make all
#	make sub

microphy_p3.o: microphy_p3.f90
routinesCld1d.o: routinesCld1d.f90
cld1d.o: cld1d.f90 routinesCld1d.o microphy_p3.o

all: $(OBJS)

obj: all

distclean: clean
	rm -f  sub

execld: execld.f90 $(OBJS)
# --> makes executable 'execld'

# Alternative:
# $ make all
# $ f90 -o execld execld.f90 *.o
# or, $ f90 -o execld execld.f90 my2mom_main_mod.o  my2mom_fncs_mod.o  my2mom_sedi_mod.o routinesCld1d.o  cld1d.o
# --> makes executable 'execld'

