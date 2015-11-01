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

.PHONY: all $(MKS) run_$(PASS) run_$(MODEL) run clean kill killf del

all: $(MKS)
	@ echo > /dev/null

$(MKS):
	@ $(MAKE) -f $(MKDIR)/$@.mk all

run: run_$(PASS)
	@ jbinfo

run_$(PASS): $(SH) $(BINDIR)/$(PASS) run_$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) && ../$< $(PROJ) $(PASS) $(MODEL) )

run_$(MODEL): $(BINDIR)/$(PASS)_$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) && ../$< $(MODELOPTS) )

run_%: $(DATDIR)/$(PASS)_%.dat | $(PWD)/$(RUNDIR)
	cp $< $(RUNDIR)/$(PASS).dat

$(PWD)/$(RUNDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINDIR) $(OBJDIR) $(DEPDIR) $(RUNDIR) $(LOGDIR)

kill killf del:
	- jbadmin -$@ -proj $(PROJ) all
