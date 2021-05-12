# Specify compiler and flags
cxx=g++-10 -fopenmp
# BLAS/LAPACK flags for linear algebra
lp_lflags=-framework Accelerate

# Lists of files to be built
objs=cir.o
src=$(patsubst %.o,%.cc,$(objs))
execs=cir_test cir_conv cir_integral

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

cir_test: cir_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

cir_conv: cir_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

cir_integral: cir_integral.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend