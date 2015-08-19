# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'model'

MAKEINC = Makefile.inc

include $(MAKEINC)

INCS =

LIBS =

LNKS = -llapacke -llapack -lcblas -lblas -ltmglib -lgfortran

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

$(BINDIR)/%: $(SRCDIR)/%.cpp $(MAKEINC) | $(PWD)/$(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCS) $(LIBS) $(LNKS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp $(MAKEINC) | $(PWD)/$(DEPDIR)
	$(CXX) -E -MM $< -MF $@ -MT '$(BINDIR)/$*' $(INCS)

$(PWD)/$(BINDIR) $(PWD)/$(DEPDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS) $(DEPS)

-include $(DEPS)
