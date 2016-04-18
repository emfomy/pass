# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile.inc

MAKEFLAGS += --no-print-directory

PROJ      ?= pass
PASS      ?= genlin
MODEL     ?= inglai
MODELOPTS ?=

CXX      = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -Wpedantic -m64

RM      = rm -rf
RANLIB  = ranlib
ARCHIVE = ar rv

# MKL
ifndef MKLROOT
  $(error MKLROOT is not set!)
endif
MKLINC   ?= -I$(MKLROOT)/include
MKLLIB   ?= -L$(MKLROOT)/lib/intel64
MKLLNK    = -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
MKLLNKSEQ = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

# MPI
ifndef MPIROOT
  $(error MPIROOT is not set!)
endif
MPIINC ?= -I$(MPIROOT)/include
MPILIB ?= -L$(MPIROOT)/lib
MPILNK  = -lmpi_cxx -lmpi

# OpenMP
OMPFLAGS = -fopenmp
