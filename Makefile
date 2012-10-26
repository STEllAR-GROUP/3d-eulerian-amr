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
        CXXFLAGS+=-ggdb -DDEBUG
    else
        CXXFLAGS+=-O3 -march=native -DNDEBUG
    endif
endif

ifdef AMR_MODULE
    CXXFLAGS+=-D$(AMR_MODULE)
endif

CXXFLAGS+=-c -I. $(OMPFLAG) 

LDFLAGS+=$(OMPFLAG) -lstdc++

ifndef SILO_ROOT
    LDFLAGS+=-L/usr/local/lib
	CXXFLAGS+=-I/usr/local/include
else
    LDFLAGS+=-L$(SILO_ROOT)/lib
    CXXFLAGS+=-I$(SILO_ROOT)/include
endif

ifdef USE_SILOH5
    LDFLAGS+=-lsiloh5
else
    LDFLAGS+=-lsilo
endif

SOURCES=$(shell find . -name '*.cpp')

OBJECTS=$(SOURCES:.cpp=.o)

SDIRS=$(shell find ! -regex '\.\(/\(build\|\.\).*\|\)' -type d -printf "%f\n")
BDIRS=$(addprefix build/, $(SDIRS))

EXECUTABLE=amr

all: directories $(SOURCES) $(EXECUTABLE)

.PHONY: directories
directories: $(BDIRS)/  

$(BDIRS)/:
	mkdir -p $@ 

$(EXECUTABLE): $(OBJECTS) 
	echo $(LDFLAGS)
	$(CC) $(addprefix build/, $(OBJECTS)) -o build/$@ $(LDFLAGS) 

.cpp.o: 
	$(CC) $(CXXFLAGS) $< -o build/$@

clean: 
	rm -rf build


