# Particle Swarm Stepwise (PaSS) Algorithm
# The Makafile for 'genlin'

MAKEINC = Makefile.inc

include $(MAKEINC)

INC =

LIB =

TARGET = genlin

BINDIR = bin

SRCDIR = src/$(TARGET)

OBJDIR = obj/$(TARGET)

BIN = $(BINDIR)/$(TARGET)

SRC = $(wildcard $(SRCDIR)/*.cpp)

HDR = $(wildcard $(SRCDIR)/*.hpp)

OBJ = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

.PHONY: all clean

all: $(BIN)
	@ echo > /dev/null

$(BIN): $(OBJ) $(MAKEINC) | $(BINDIR) 
	$(ARCHIVE) $@ $(OBJ)
	- $(RANLIB) $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HDR) $(MAKEINC) | $(OBJDIR) 
	$(BGCXX) $(BGCXXFLAGS) -c $< $(INC) $(INCLUDE) -o $@

$(BINDIR) $(OBJDIR):
	mkdir -p $@

clean:
	$(RM) $(BIN) $(OBJ)
