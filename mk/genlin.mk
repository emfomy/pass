# Particle Swarm Stepwise (PaSS) Algorithm
# The Makefile for 'genlin'

include Makefile.inc

CXXFLAGS += $(OMPFLAGS)
INCS = $(MPIINC) $(MKLINC) $(MAGMAINC) $(CUDAINC)
LIBS = $(MPILIB) $(MKLLIB) $(MAGMALIB) $(CUDALIB)
LNKS = $(MPILNK) $(MKLLNK) $(MAGMALNK) $(CUDALNK)

NAME = genlin
BINDIR = bin
SRCDIR = src/$(NAME)
OBJDIR = obj/$(NAME)
DEPDIR = dep/$(NAME)

BINS = $(BINDIR)/$(NAME)
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
DEPS = $(OBJS:$(OBJDIR)/%.o=$(DEPDIR)/%.d)

.PHONY: all dep clean

all: $(BINS)
	@ echo > /dev/null

dep: $(DEPS)
	@ echo > /dev/null

$(BINS): $(OBJS) | $(PWD)/$(BINDIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LIBS) $(LNKS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(PWD)/$(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp | $(PWD)/$(DEPDIR)
	@ $(CXX) $(CXXFLAGS) -E -MM $< -MF $@ -MT '$(OBJDIR)/$*.o' $(INCS)

$(PWD)/$(BINDIR) $(PWD)/$(OBJDIR) $(PWD)/$(DEPDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS) $(OBJS) $(DEPS)

-include $(DEPS)
