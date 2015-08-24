# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile.inc

MAKEFLAGS += --no-print-directory

PROJ      = pass
PASS      = genlin
MODEL     = inglai
MODELOPTS =

CXX      = g++
CXXFLAGS = -O3 -std=c++11

RM      = rm -rf
RANLIB  = ranlib
ARCHIVE = ar rv

MPILNK = -lmpi_cxx -lmpi
MKLLNK = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
OMPLNK = -fopenmp
