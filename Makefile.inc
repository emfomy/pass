# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile.inc

MAKEFLAGS += --no-print-directory

PROJ      ?= pass
PASS      ?= genlin
MODEL     ?= inglai
PASSOPTS  ?=
MODELOPTS ?=

CXX      = mpicxx
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -Wpedantic -m64

RM      = rm -rf
RANLIB  = ranlib
ARCHIVE = ar rv

# OpenMP
OMPFLAGS = -fopenmp

# MPI
ifndef MPIROOT
  $(error MPIROOT is not set!)
endif
MPIINC ?= -I$(MPIROOT)/include
MPILIB ?= -L$(MPIROOT)/lib
MPILNK  = -lmpi_cxx -lmpi

# CUDA
ifndef CUDAROOT
  $(error CUDAROOT is not set!)
endif
CUDAINC ?= -isystem $(CUDAROOT)/include
CUDALIB ?= -L$(CUDAROOT)/lib64
CUDALNK  = -lcudart -lcublas -lcusparse

# MKL
ifndef MKLROOT
  $(error MKLROOT is not set!)
endif
MKLINC ?= -I$(MKLROOT)/include
MKLLIB ?= -L$(MKLROOT)/lib/intel64
MKLLNK  = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

# MAGMA
ifndef MAGMAROOT
  $(error MAGMAROOT is not set!)
endif
MAGMAINC ?= -I$(MAGMAROOT)/include
MAGMALIB ?= -L$(MAGMAROOT)/lib
MAGMALNK  = -lmagma
