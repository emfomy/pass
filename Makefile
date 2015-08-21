# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile

MAKEINC = Makefile.inc

include $(MAKEINC)

PROJ      = pass
PASS      = genlin
MODEL     = power
MODELOPTS =

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

run: .$(PASS)
	@ echo > /dev/null

.%: $(SH) $(BINDIR)/% .$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< $(PROJ) $(PASS) $(MODEL) )

.%: $(BINDIR)/$(PASS)_% | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< $(MODELOPTS) )

.%: $(DATDIR)/$(PASS)_%.dat | $(PWD)/$(RUNDIR)
	cp $< $(RUNDIR)/$(PASS).dat

$(PWD)/$(RUNDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINDIR) $(OBJDIR) $(DEPDIR) $(RUNDIR) $(LOGDIR)

kill killf del:
	- jbadmin -$@ -proj $(PROJ) all
