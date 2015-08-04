# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile.inc

MAKEFLAGS += --no-print-directory

CXX        = g++
CXXFLAGS   = -O3 -std=c++0x
BGCXX      = mpixlcxx_r
BGCXXFLAGS = -O3 -qlanglvl=extended0x -qsmp=omp -qstrict

RM      = rm -rf
RANLIB  = ranlib
ARCHIVE = ar rv

INCLUDE = \
	-I$(ESSLINC)

LIBRARY = \
	-L$(ESSLLIB) -lesslbg \
	-L$(XLFLIB) -lxlopt -lxlf90_r -lxlfmath -lxl
