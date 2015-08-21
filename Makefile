# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile

MAKEINC = Makefile.inc

include $(MAKEINC)

PROJ  = pass
NAME  = genlin
MODEL = hung_1_1

SRCDIR = src
BINDIR = bin
OBJDIR = obj
DEPDIR = dep
RUNDIR = run
LOGDIR = log
DATDIR = dat
MKDIR  = mk
SHDIR  = sh

SH = $(SHDIR)/pass.sh

MKS = $(notdir $(basename $(wildcard $(MKDIR)/*.mk)))

.PHONY: all $(MKS) run clean cancel

all: $(MKS)
	@ echo > /dev/null

$(MKS):
	@ $(MAKE) -f $(MKDIR)/$@.mk all

run: .$(NAME)
	@ echo > /dev/null

.%: $(SH) $(BINDIR)/% .$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< $(PROJ) $* )

.%: $(BINDIR)/$(NAME)_% | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< )

.%: $(DATDIR)/$(NAME)_%.dat | $(PWD)/$(RUNDIR)
	cp $< $(RUNDIR)/$(NAME).dat

$(PWD)/$(RUNDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINDIR) $(OBJDIR) $(DEPDIR) $(RUNDIR) $(LOGDIR)

kill killf del:
	- jbadmin -$@ -proj $(PROJ) all
