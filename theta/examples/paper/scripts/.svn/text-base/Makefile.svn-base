CXX := g++

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTINCDIR   := $(shell root-config --incdir)

CXXFLAGS := -g -Wall -fPIC $(ROOTCFLAGS)

LIBS := $(ROOTLIBS) -lm

all: make_histos

%: %.cc
	@echo -n "Compiling $< ..."
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@
	@echo "Done"

