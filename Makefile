# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile

MAKEINC = Makefile.inc

include $(MAKEINC)

TGTDIR = mk

RUNDIR = run

MKS = $(notdir $(basename $(wildcard $(TGTDIR)/*.mk)))

MODEL = inglai

DEMO = sh/genlin.sh

.PHONY: all $(MKS) run clean cancel

all: $(MKS)
	@ echo > /dev/null

$(MKS):
	@ $(MAKE) -f $(TGTDIR)/$@.mk all

run: $(DEMO) .$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< )

.%: bin/genlin_% | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) ; ../$< )

.%: data/genlin_%.dat | $(PWD)/$(RUNDIR)
	cp $< $(RUNDIR)/genlin.dat

$(PWD)/$(RUNDIR):
	@ mkdir -p $@

clean:
	$(RM) bin obj dep run

cancel:
	scancel -A $(USER)
