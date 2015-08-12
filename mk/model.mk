# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'model'

MAKEINC = Makefile.inc

include $(MAKEINC)

INCS =

LIBS = -llapacke -llapack -lcblas -lblas -ltmglib -lgfortran

TGTDIR = model

BINDIR = bin

SRCDIR = src/$(TGTDIR)

DEPDIR = dep/$(TGTDIR)

SRCS = $(wildcard $(SRCDIR)/*.cpp)

BINS = $(SRCS:$(SRCDIR)/%.cpp=$(BINDIR)/%)

DEPS = $(BINS:$(BINDIR)/%=$(DEPDIR)/%.d)

.PHONY: all dep run clean

all: $(BINS)
	@ echo > /dev/null

dep: $(DEPS)
	@ echo > /dev/null

$(BINDIR)/%: $(SRCDIR)/%.cpp $(MAKEINC) | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCS) $(LIBS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp $(MAKEINC) | $(DEPDIR)
	$(CXX) -E -MM $< -MF $@ -MT '$(BINDIR)/$*' $(INCS)

$(BINDIR) $(DEPDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS) $(DEPS)

-include $(DEPS)
