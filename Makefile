# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile

MAKEINC = Makefile.inc

include $(MAKEINC)

TGTDIR = mk

TGT = $(wildcard $(TGTDIR)/*.mk)

MODEL = src/data/genlin_hung_1_1.dat

RUN = sh/genlin.sh

.PHONY: all $(TGT) run clean cancel

all: $(TGT)
	@ echo > /dev/null

$(TGT):
	@ ( $(MAKE) -f $@ all )

run: $(MODEL) $(RUN)
	@ mkdir -p $@
	cp $(MODEL) run/genlin.dat
	( cd run ; ../$(RUN) )

clean:
	@ for tgt in $(TGT) ; do ( $(MAKE) -f $$tgt clean ) done
	$(RM) bin obj run

cancel:
	scancel -A xntmuyang
