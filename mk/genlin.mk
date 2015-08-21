# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'genlin'

MAKEINC = Makefile.inc

include $(MAKEINC)

INCS = $(MKLINC) $(MPIINC)
LIBS = $(MKLLIB) $(MPILIB)
LNKS = $(MKLLNK) $(MPILNK) $(OMPLNK)

NAME = genlin
BINDIR = bin
SRCDIR = src/$(NAME)
OBJDIR = obj/$(NAME)
DEPDIR = dep/$(NAME)

BINS = $(BINDIR)/$(NAME)
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
DEPS = $(OBJS:$(OBJDIR)/%.o=$(DEPDIR)/%.d)

.PHONY: all dep run clean

all: $(BINS)
	@ echo > /dev/null

dep: $(DEPS)
	@ echo > /dev/null

$(BINS): $(OBJS) $(MAKEINC) | $(PWD)/$(BINDIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LIBS) $(LNKS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(MAKEINC) | $(PWD)/$(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp $(MAKEINC) | $(PWD)/$(DEPDIR)
	@ $(CXX) $(CXXFLAGS) -E -MM $< -MF $@ -MT '$(OBJDIR)/$*.o' $(INCS)

$(PWD)/$(BINDIR) $(PWD)/$(OBJDIR) $(PWD)/$(DEPDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS) $(OBJS) $(DEPS)

-include $(DEPS)
