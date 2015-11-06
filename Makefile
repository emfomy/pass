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
HTMLDIR = html
PASSHTML = pass.html

SH = $(SHDIR)/pass.sh

MKS = $(notdir $(basename $(wildcard $(MKDIR)/*.mk)))
DOCS = $(notdir $(basename $(wildcard $(DOCDIR)/*.inc)))

.PHONY: all $(MKS) doc $(DOCS) run run-$(PASS) run-$(MODEL) clean kill killf del

all: $(MKS)
	@ echo > /dev/null

$(MKS):
	@ $(MAKE) -f $(MKDIR)/$@.mk all

doc: $(DOCS)
	@ sed -i '' 's|PaSS Documentation|Particle Swarm Stepwise (PaSS) Algorithm|g' $(HTMLDIR)/index.html
	ln -sf $(HTMLDIR)/index.html $(PASSHTML)

$(DOCS): | $(PWD)/$(HTMLDIR)
	doxygen $(DOCDIR)/$@.inc

run: run-$(PASS)
	@ jbinfo

run-$(PASS): $(SH) $(BINDIR)/$(PASS) run-$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) && ../$< $(PROJ) $(PASS) $(MODEL) )

run-$(MODEL): $(BINDIR)/$(PASS)_$(MODEL) | $(PWD)/$(RUNDIR)
	( cd $(RUNDIR) && ../$< $(MODELOPTS) )

$(PWD)/$(RUNDIR) $(PWD)/$(HTMLDIR):
	@ mkdir -p $@

clean:
	$(RM) $(BINDIR) $(OBJDIR) $(DEPDIR) $(HTMLDIR) $(PASSHTML) $(RUNDIR) $(LOGDIR)

kill killf del:
	- jbadmin -$@ -proj $(PROJ) all
