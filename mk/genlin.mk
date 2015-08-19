# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'genlin'

MAKEINC = Makefile.inc

include $(MAKEINC)

INCS = -I$(ESSLINC)

LIBS = -L$(ESSLLIB) -L$(XLFLIB)

LNKS = -lesslbg -lxlopt -lxlf90_r -lxlfmath -lxl

TGTDIR = genlin

BINDIR = bin

SRCDIR = src/$(TGTDIR)

OBJDIR = obj/$(TGTDIR)

DEPDIR = dep/$(TGTDIR)

BINS = $(BINDIR)/$(TGTDIR)

SRCS = $(wildcard $(SRCDIR)/*.cpp)

OBJS = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

DEPS = $(OBJS:$(OBJDIR)/%.o=$(DEPDIR)/%.d)

.PHONY: all dep run clean

all: $(BINS)
	@ echo > /dev/null

dep: $(DEPS)
	@ echo > /dev/null

$(BINS): $(OBJS) $(MAKEINC) | $(PWD)/$(BINDIR)
	$(BGCXX) $(BGCXXFLAGS) $(OBJS) -o $@ $(LIBS) $(LNKS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(MAKEINC) | $(PWD)/$(OBJDIR)
	$(BGCXX) $(BGCXXFLAGS) -c $< -o $@ $(INCS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp $(MAKEINC) | $(PWD)/$(DEPDIR)
	$(CXX) -E -MM $< -MF $@ -MT '$(OBJDIR)/$*.o' $(INCS)

$(PWD)/$(BINDIR) $(PWD)/$(OBJDIR) $(PWD)/$(DEPDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS) $(OBJS) $(DEPS)

-include $(DEPS)
