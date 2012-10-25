.SUFFIXES: .cpp .c .f90
CP=icpc
CC=icc
FC=ifort
LD=icpc
FLAGS=-I. -I/home/dmarce1/include -DNDEBUG -fast -openmp 
#FLAGS=-I. -I/home/dmarce1/include -DNDEBUG -pg -openmp-stubs -O0
#FLAGS=-I.  -DNDEBUG -g -O0 -openmp
LDLIBS=-lsilo -L/home/dmarce1/lib 
CPPSOURCES=assert.cpp indexer_2d.cpp program.cpp indexer_2d_by2.cpp main.cpp reconstruct.cpp grid/grid_amr_bound.cpp grid/grid_interface.cpp grid/grid_phys_bound.cpp grid/grid.cpp grid/grid_node.cpp binary/binary.cpp binary/initialize.cpp binary/binary_driver.cpp  binary/state.cpp binary/state_ztwd.cpp poisson/poisson_amr_bound.cpp poisson/poisson_phys_bound.cpp poisson/poisson.cpp euler/hydro_driver.cpp euler/initialize.cpp euler/state.cpp oct_node/child_index.cpp oct_node/oct_face.cpp oct_node/oct_node.cpp helmholtz/helmeos.cpp
CSOURCES=binary/lane_emden.c 

CPPOBJECTS=$(CPPSOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.c=.o)
EXECUTABLE=amr

.cpp.o:
	$(CP) $(FLAGS) -c $< -o $@

.c.o:
	$(CC) $(FLAGS) -c $< -o $@


$(EXECUTABLE): $(CPPOBJECTS) $(COBJECTS) 
	$(LD) $(FLAGS) $(CPPOBJECTS) $(COBJECTS) -o $@ $(LDLIBS)

all: $(CPPSOURCES) $(FSOURCES) $(CSOURCES) $(EXECUTABLE)
	
clean:
	rm *.o
	rm ./*/*.o
	rm amr

