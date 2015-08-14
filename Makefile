# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile

MAKEINC = Makefile.inc

include $(MAKEINC)

TGTDIR = mk

RUNDIR = run

TGTS = $(notdir $(basename $(wildcard $(TGTDIR)/*.mk)))

MODEL = inglai

MAIN = sh/genlin.sh

.PHONY: all $(TGTS) demo clean cancel

all: $(TGTS)
	@ echo > /dev/null

$(TGTS):
	@ ( $(MAKE) -f $(TGTDIR)/$@.mk dep all )

demo: $(MAIN) .$(MODEL) | $(RUNDIR)
	( cd $(RUNDIR) ; ../$< )

.%: bin/genlin_% | $(RUNDIR)
	( cd $(RUNDIR) ; ../$< )

.%: data/genlin_%.dat | $(RUNDIR)
	cp $< $(RUNDIR)/genlin.dat

$(RUNDIR):
	@ mkdir -p $@

clean:
	$(RM) bin obj dep run

cancel:
	scancel -A $(USER)
