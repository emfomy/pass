# Particle Swarm Stepwise (PaSS) Algorithm
# The Makefile for 'data'

include Makefile.inc

BINDIR = bin
DATDIR = dat

DATS = $(wildcard $(DATDIR)/*.dat)
BINS = $(DATS:$(DATDIR)/%.dat=$(BINDIR)/%)

.PHONY: all clean

all: $(BINS)
	@ echo > /dev/null

$(BINDIR)/%: $(DATDIR)/%.dat | $(PWD)/$(BINDIR)
	printf "#!/bin/bash\n\ncp $(PWD)/$< $(PASS).dat\n" > $@ ; chmod +x $@

$(PWD)/$(BINDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINS)
