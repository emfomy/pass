# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'model'

MAKEINC = Makefile.inc

include $(MAKEINC)

INC = -I$(CBLASINC) -I$(LAPACKINC)

LIB = -L$(CBLASLIB) -L$(LAPACKLIB) -lcblas -llapacke -llapack -lrefblas -ltmglib -lgfortran

TARGET = model

BINDIR = bin

SRCDIR = src/$(TARGET)

SRC = $(wildcard $(SRCDIR)/*.cpp)

BIN = $(SRC:$(SRCDIR)/%.cpp=$(BINDIR)/%)

.PHONY: all run clean

all: $(BIN)
	@ echo > /dev/null

$(BINDIR)/%: $(SRCDIR)/%.cpp $(MAKEINC) | $(BINDIR) 
	$(CXX) $(CXXFLAGS) $< $(INC) $(LIB) -o $@

$(BINDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BIN)
