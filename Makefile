
MAKEFILE_OPT = $(PWD)/Makefile_opt.in
include $(MAKEFILE_OPT)

MAKEFILE_H5  = $(PWD)/Makefile_dir.in
include $(MAKEFILE_H5)

TEMPLATES = bexp bx3d earth fieldloop flock flock_grmhd isentropic jupiter kepler kh mri2 rotor shear shocktube spinring spread vortex sorathia_grmhd blast_grmhd fieldloop_grmhd bl kep_ring torus_fm isentropic_RAM entropywave acousticwave cb acousticwave_cart bondi binarybondi alfvenwave alfvenwave_cart magnetosonicwave_cart magnetosonicwave_cart3d bondi_mhd fieldloop_linear_cart ecc tde acousticwave_sph cartesian_shear cartesian_shear_sph vortexShedding cb_cc cb_cc3d cb_amv singleSinkTest

GIT_VERSION = $(shell git describe --dirty --always --tags)

# Disable cartesian interpolation by default.
ifndef ENABLE_CART_INTERP
ENABLE_CART_INTERP = 0
endif

GIT_VERSION := $(strip $(GIT_VERSION))
INITIAL := $(strip $(INITIAL))
HYDRO := $(strip $(HYDRO))
GEOMETRY := $(strip $(GEOMETRY))
BOUNDARY := $(strip $(BOUNDARY))
OUTPUT := $(strip $(OUTPUT))
RESTART := $(strip $(RESTART))
PLANETS := $(strip $(PLANETS))
HLLD := $(strip $(HLLD))
ANALYSIS := $(strip $(ANALYSIS))
REPORT := $(strip $(REPORT))
METRIC := $(strip $(METRIC))
FRAME := $(strip $(FRAME))
ENABLE_CART_INTERP := $(strip $(ENABLE_CART_INTERP))
NUM_C := $(strip $(NUM_C))
NUM_N := $(strip $(NUM_N))
CT_MODE := $(strip $(CT_MODE))
USE_MPI := $(strip $(USE_MPI))

OPT_DEFS = -DGIT_VERSION=\"$(GIT_VERSION)\"
OPT_DEFS += -DINITIAL=\"$(INITIAL)\"
OPT_DEFS += -DHYDRO=\"$(HYDRO)\"
OPT_DEFS += -DGEOMETRY=\"$(GEOMETRY)\"
OPT_DEFS += -DBOUNDARY=\"$(BOUNDARY)\"
OPT_DEFS += -DOUTPUT=\"$(OUTPUT)\"
OPT_DEFS += -DRESTART=\"$(RESTART)\"
OPT_DEFS += -DPLANETS=\"$(PLANETS)\"
OPT_DEFS += -DHLLD=\"$(HLLD)\"
OPT_DEFS += -DANALYSIS=\"$(ANALYSIS)\"
OPT_DEFS += -DREPORT=\"$(REPORT)\"
OPT_DEFS += -DMETRIC=\"$(METRIC)\"
OPT_DEFS += -DFRAME=\"$(FRAME)\"
OPT_DEFS += -DENABLE_CART_INTERP=$(ENABLE_CART_INTERP)
OPT_DEFS += -DNUM_C=$(NUM_C)
OPT_DEFS += -DNUM_N=$(NUM_N)
OPT_DEFS += -DCT_MODE=$(CT_MODE)

DIR_DEFS = -DUSE_MPI=$(USE_MPI)

DEBUG_FLAGS = -g
WARNING_FLAGS = -Wall #-Wextra -Wshadow -Wconversion -Werror

FLAGS = -O3 $(OPT_DEFS) $(DIR_DEFS) $(DEBUG_FLAGS) $(WARNING_FLAGS)

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lhdf5 -lm

#INC = -I/usr/include
#LIB = -L/usr/include -lhdf5 -lm

OBJ = main.o readpar.o timestep.o onestep.o riemann.o mpisetup.o gridsetup.o domain.o misc.o $(GEOMETRY).o faces_alt.o exchange.o plm.o report.o profiler.o planet.o omega.o analysis.o $(HLLD).o rotframe.o boundary_functions.o geometry_functions.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o $(BOUNDARY).o $(RESTART).o $(PLANETS).o $(METRIC).o $(FRAME).o calc.a $(ANALYSIS).o $(REPORT).o noise.o sink.o #snapshot.o

CALC_OBJ = Calc/bondi.o Calc/integrate.o Calc/magnetosonic.o
CT_OBJ = bfields.o

ifneq ($(CT_MODE), 0)
OBJ += $(CT_OBJ)
endif


default: disco

.PHONY: $(TEMPLATES)

$(TEMPLATES):
	cp Templates/$@.par in.par
	cp Templates/$@.in Makefile_opt.in
	make clean
	make

%.o: %.c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c $< -o $@

calc.a: $(CALC_OBJ)
	ar rcs $@ $^

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Hydro/$(HYDRO).c

$(GEOMETRY).o : Geometry/$(GEOMETRY).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Geometry/$(GEOMETRY).c

$(PLANETS).o : Planets/$(PLANETS).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Planets/$(PLANETS).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Output/$(OUTPUT).c

$(RESTART).o : Restart/$(RESTART).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Restart/$(RESTART).c

$(ANALYSIS).o : Diagnostics/$(ANALYSIS).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Diagnostics/$(ANALYSIS).c

$(REPORT).o : Reports/$(REPORT).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Reports/$(REPORT).c

$(METRIC).o : Hydro/Metric/$(METRIC).c paul.h Hydro/metric.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Hydro/Metric/$(METRIC).c

$(FRAME).o : Hydro/Frame/$(FRAME).c paul.h Hydro/frame.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Hydro/Frame/$(FRAME).c

disco: $(OBJ) paul.h
	$(CC) $(FLAGS) -o disco $(OBJ) $(LOCAL_LDFLAGS) $(LIB)

clean:
	rm -f Calc/*.o *.a *.o disco
