
CC=icc
CFLAGS=-c -I. -openmp -fast -D$(AMR_MODULE) -DNDEBUG
LDFLAGS=-openmp -L/usr/local/lib -fast -lsilo -lstdc++
SOURCES=assert.cpp indexer_2d.cpp program.cpp indexer_2d_by2.cpp main.cpp reconstruct.cpp grid/grid_amr_bound.cpp grid/grid_interface.cpp grid/grid_phys_bound.cpp grid/grid.cpp grid/grid_node.cpp binary/binary.cpp binary/initialize.cpp binary/binary_driver.cpp binary/lane_emden.c binary/state.cpp poisson/poisson_amr_bound.cpp poisson/poisson_phys_bound.cpp poisson/poisson.cpp euler/hydro_driver.cpp euler/initialize.cpp euler/state.cpp oct_node/child_index.cpp oct_node/oct_face.cpp oct_node/oct_node.cpp 


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=amr

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm *.o
	rm ./*/*.o
	rm amr


