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
DOCDIR = doc
HTMLDIR = docs
PASSHTML = pass.html

SH = $(SHDIR)/pass.sh

MKS = $(notdir $(basename $(wildcard $(MKDIR)/*.mk)))
DOCS = $(notdir $(basename $(wildcard $(DOCDIR)/*.inc)))

.PHONY: all $(MKS) doc $(DOCS) run run-pass run-model run-$(PASS) run-$(MODEL) clean kill killf del

all: $(MKS)
	@ echo > /dev/null

$(MKS):
	@ $(MAKE) -f $(MKDIR)/$@.mk all

doc: $(DOCS) $(PASSHTML)
	sed -ig 's|PaSS Documentation|Particle Swarm Stepwise (PaSS) Algorithm|g' $(HTMLDIR)/index.html

$(DOCS): | $(PWD)/$(HTMLDIR)
	doxygen $(DOCDIR)/$@.inc

$(PASSHTML):
	echo "<html><META HTTP-EQUIV='refresh' CONTENT='0; URL=$(HTMLDIR)/index.html'></html>" > $@

run: $(MKS) run-$(PASS)
	@ jbinfo

run-pass: run-$(PASS)

run-model: run-$(MODEL)

run-$(PASS): $(SH) $(BINDIR)/$(PASS) run-$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) && ../$< $(PROJ) $(PASS) $(MODEL) $(PASSOPTS) )

run-$(MODEL): $(BINDIR)/$(PASS)_$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) && ../$< $(MODELOPTS) )

$(PWD)/$(RUNDIR) $(PWD)/$(HTMLDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINDIR) $(OBJDIR) $(DEPDIR) $(RUNDIR)

distclean:
	$(RM) $(BINDIR) $(OBJDIR) $(DEPDIR) $(HTMLDIR) $(PASSHTML) $(RUNDIR) $(LOGDIR)

kill killf del:
	- jbadmin -$@ -proj $(PROJ) all
