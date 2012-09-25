ifdef CC
    ifneq ($(origin CC),'command line') 
        CC=icc 
    endif
endif

ifeq ($(CC),icc)
    OMPFLAG=-openmp
    ifdef DEBUG
        CXXFLAGS+=-g -DDEBUG
    else
        CXXFLAGS+=-fast -DNDEBUG
    endif
else 
    OMPFLAG=-fopenmp
    ifdef DEBUG
        CXXFLAGS+=-O0 -ggdb -DDEBUG -fno-inline
    else
        CXXFLAGS+=-O3 -march=native -DNDEBUG
    endif
endif

ifdef AMR_MODULE
    CXXFLAGS+=-D$(AMR_MODULE)
endif

CXXFLAGS=-c -I. $(OMPFLAG) 
LDFLAGS=$(OMPFLAG) -L/usr/local/lib -lsiloh5 -lstdc++

SOURCES=assert.cpp \
        indexer_2d.cpp \
        program.cpp \
        indexer_2d_by2.cpp \
        main.cpp \
        reconstruct.cpp \
        grid/grid_amr_bound.cpp \
        grid/grid_interface.cpp \
        grid/grid_phys_bound.cpp \
        grid/grid.cpp \
        grid/grid_node.cpp \
        binary/binary.cpp \
        binary/initialize.cpp \
        binary/binary_driver.cpp \
        binary/lane_emden.cpp \
        binary/state.cpp \
        poisson/poisson_amr_bound.cpp \
        poisson/poisson_phys_bound.cpp \
        poisson/poisson.cpp \
        euler/hydro_driver.cpp \
        euler/initialize.cpp \
        euler/state.cpp \
        oct_node/child_index.cpp \
        oct_node/oct_face.cpp \
        oct_node/oct_node.cpp

OBJECTS=$(SOURCES:.cpp=.o)

DIRECTORIES=build \
            build/grid \
            build/binary \
            build/poisson \
            build/euler \
            build/oct_node

EXECUTABLE=amr

all: directories $(SOURCES) $(EXECUTABLE)

.PHONY: directories
directories: $(DIRECTORIES)/  

$(DIRECTORIES)/:
	mkdir -p $@ 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(addprefix build/, $(OBJECTS)) -o build/$@ $(LDFLAGS) 

.cpp.o: 
	$(CC) $(CXXFLAGS) $< -o build/$@

clean: 
	rm -rf build


