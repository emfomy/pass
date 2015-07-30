# Particle Swarm Stepwise (PaSS) Algorithm
# The main Makefile

TOP = .

MAKEINC = $(TOP)/Makefile.inc

include $(MAKEINC)

TGTDIR = mk

TGT = $(wildcard $(TGTDIR)/*.mk)

.PHONY: all $(TGT) run clean

all: $(TGT)
	@ echo > /dev/null

$(TGT):
	@ ( $(MAKE) -f $@ all )

clean:
	@ for tgt in $(TGT) ; do ( $(MAKE) -f $$tgt clean ) done
	$(RM) bin obj run
