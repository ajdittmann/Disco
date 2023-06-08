MAKEFILE_H5 = $(PWD)/../Makefile_dir.in
include $(MAKEFILE_H5)

EXE  = vdisco
RNDR = rdisco
ALL = $(EXE) $(RNDR)

GLUT = /usr/X11R6
HDF5 = $(H55)

CC_SDL := $(shell sdl2-config --libs)
LD_SDL := $(shell sdl2-config --cflags)

LIBS = $(LOCAL_LDFLAGS)
CFLAGS = $(LOCAL_CFLAGS) -Wall -Wextra -g -O3

#LIBS = -L$(GLUT)/lib -L$(HDF5)/lib -framework GLUT -framework OpenGL -lXmu -lXext -lXi -lX11 -lm -lhdf5 -lz
UNAME = $(shell uname -s)
ifeq ($(UNAME), Linux)
	GLUT = /usr/include
	LIBS += -L$(HDF5)/lib -L$(GLUT)/lib -lGL -lGLU -lglut -lhdf5 -lm
	CFLAGS += -D LINUX
endif
ifeq ($(UNAME), Darwin)
	LIBS += -L$(HDF5)/lib -framework GLUT -framework OpenGL -lhdf5 -lm
	CFLAGS += -D OSX
endif

INCL = -I$(GLUT)/include -I$(HDF5)/include
LIBS += -lhdf5 -lm

CC = gcc

default: $(ALL)

$(EXE): $(EXE).c colorbar.c gifer.c gifer.h vutil.c vutil.h
	$(CC) -c colorbar.c
	$(CC) -c gifer.c
	$(CC) -c vutil.c $(INCL)
	$(CC) $(CFLAGS) -o $@ $< colorbar.o gifer.o vutil.o $(INCL) $(LIBS)

$(RNDR): $(RNDR)_sdl.c rdisco_core.c colorbar.c gifer.c gifer.h vutil.c vutil.h rdisco.h
	$(CC) -c vutil.c $(INCL)
	$(CC) -c rdisco_core.c $(INCL)
	$(CC) $(CFLAGS) $(CC_SDL) -o $@ $< vutil.o rdisco_core.o $(INCL) $(LD_SDL) $(LIBS)

clean:
	rm -f *.o $(EXE) $(RNDR)
