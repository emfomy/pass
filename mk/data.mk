# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'data'

MAKEINC = Makefile.inc

include $(MAKEINC)

BINDIR = bin
DATDIR = dat

DATS = $(wildcard $(DATDIR)/*.dat)
BINS = $(DATS:$(DATDIR)/%.dat=$(BINDIR)/%)

.PHONY: all clean

all: $(BINS)
	@ echo > /dev/null

$(BINDIR)/%: $(DATDIR)/%.dat $(MAKEINC) | $(PWD)/$(BINDIR)
	echo -e "#!/bin/bash\n\ncp $(PWD)/$< $(PASS).dat\n" > $@
	chmod +x $@

$(PWD)/$(BINDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS)
