# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile.inc

MAKEFLAGS += --no-print-directory

PROJ      ?= pass
PASS      ?= genlog
MODEL     ?= inglai
MODELOPTS ?=

CXX      = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -pedantic

RM      = rm -rf
RANLIB  = ranlib
ARCHIVE = ar rv

# MKL
ifndef MKLROOT
  $(error MKLROOT is not set!)
endif
MKLINC ?= -I$(MKLROOT)/include -I$HOME/local/include/intel/mkl
MKLLIB ?= -L$(MKLROOT)/lib/intel64
MKLLNK  = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

# MPI
ifndef MPIROOT
  $(error MPIROOT is not set!)
endif
MKLINC ?= -I$(MPIROOT)/include
MKLLIB ?= -L$(MPIROOT)/lib
MPILNK  = -lmpi_cxx -lmpi

# OpenMP
OMPFLAGS = -fopenmp
