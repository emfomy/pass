# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile

MAKEINC = Makefile.inc

include $(MAKEINC)

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

.PHONY: all $(MKS) .$(PASS) .$(MODEL) run clean cancel

all: $(MKS)
	@ echo > /dev/null

$(MKS):
	@ $(MAKE) -f $(MKDIR)/$@.mk all

run: .$(PASS)
	@ jbinfo

.$(PASS): $(SH) $(BINDIR)/$(PASS) .$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< $(PROJ) $(PASS) $(MODEL) )

.$(MODEL): $(BINDIR)/$(PASS)_$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< $(MODELOPTS) )

.%: $(DATDIR)/$(PASS)_%.dat | $(PWD)/$(RUNDIR)
	cp $< $(RUNDIR)/$(PASS).dat

$(PWD)/$(RUNDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINDIR) $(OBJDIR) $(DEPDIR) $(RUNDIR) $(LOGDIR)

kill killf del:
	- jbadmin -$@ -proj $(PROJ) all
