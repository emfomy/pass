# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'genlin'

MAKEINC = Makefile.inc

include $(MAKEINC)

INCS = \
	-I$(ESSLINC)

LIBS = \
	-L$(ESSLLIB) -lesslbg \
	-L$(XLFLIB) -lxlopt -lxlf90_r -lxlfmath -lxl

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

$(BINS): $(OBJS) $(MAKEINC) | $(BINDIR)
	$(BGCXX) $(BGCXXFLAGS) $(OBJS) -o $@ $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(MAKEINC) | $(OBJDIR)
	$(BGCXX) $(BGCXXFLAGS) -c $< -o $@ $(INCS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp $(MAKEINC) | $(DEPDIR)
	$(CXX) -E -MM $< -MF $@ -MT '$(OBJDIR)/$*.o' $(INCS)

$(BINDIR) $(OBJDIR) $(DEPDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS) $(OBJS) $(DEPS)

-include $(DEPS)
