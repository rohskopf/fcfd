.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for the program GEO.
# Please modify the variables properly.
# We recommend to use gcc >= 4.8 with openMPI >= 1.8
#-----------------------------------------------

# Enter the path to your LAMMPS directory
LAMMPSDIR = /home/drew/lammps-5Jun19/src
#LAMMPSDIR = /home/rohskopf/lammps/src

# Enter path to the shared library (.so) file
LAMMPSLIB = /home/drew/lammps-5Jun19/src/liblammps_g++_openmpi.so
#LAMMPSLIB = /home/rohskopf/lammps/src/liblammps_g++_openmpi.so

#-----------------------------------------------

CXX = mpicxx
CXXFLAGS = -std=c++11 -Wno-write-strings

INCLUDE = -I${LAMMPSDIR}

CXXL =
LDFLAGS = 

LAPACK =
LIBS =

#-----------------------------------------------
# General rules 
#-----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = fcfd

CXXSRC= geo.cpp main.cpp mem.cpp in.cpp verify.cpp ifc2mcc.cpp asr.cpp compute.cpp

OBJS= ${CXXSRC:.cpp=.o}

default: fcfd

all: ${PROG}

fcfd: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${LIBS} ${LAMMPSLIB}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
