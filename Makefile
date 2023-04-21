PROG =	riemann_rmhd

SRCS =	Interfaces.f90 initialdata.f90 quartic.f90\
nrutil.f90 lubksb.f90 ludcmp.f90 postshock.f90 riemann_rmhd.f90

OBJS =	Interfaces.o initialdata.o quartic.o\
nrutil.o lubksb.o ludcmp.o postshock.o riemann_rmhd.o

LIBS =	

F90 = gfortran -ffree-line-length-none

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(OBJ_TEST) *.kmo *.mod

flush:
	rm -f $(PROG) $(OBJS) $(OBJ_TEST) $(OBJ_TESTV) *.kmo *.mod

data:	
	rm -f *.dat *.sol

tar:
	tar -czvf rmhd_riemann_solver.tgz $(SRCS) Makefile README RInput.txt

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

# DO NOT DELETE
