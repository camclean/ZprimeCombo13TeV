ARCH_LOC_1 := $(wildcard $(shell root-config --prefix)/test/Makefile.arch)
ARCH_LOC_2 := $(wildcard $(shell root-config --prefix)/share/root/test/Makefile.arch)
ARCH_LOC_3 := $(wildcard $(shell root-config --prefix)/share/doc/root/test/Makefile.arch)
ARCH_LOC_4 := $(wildcard $(shell root-config --prefix)/etc/Makefile.arch)
ARCH_LOC_5 := $(wildcard $(shell root-config --prefix)/etc/root/Makefile.arch)
ifneq ($(strip $(ARCH_LOC_1)),)
  $(info Using $(ARCH_LOC_1))
  include $(ARCH_LOC_1)
else
  ifneq ($(strip $(ARCH_LOC_2)),)
    $(info Using $(ARCH_LOC_2))
    include $(ARCH_LOC_2)
  else
    ifneq ($(strip $(ARCH_LOC_3)),)
      $(info Using $(ARCH_LOC_3))
      include $(ARCH_LOC_3)
    else
                ifneq ($(strip $(ARCH_LOC_4)),)
        $(info Using $(ARCH_LOC_4))
        include $(ARCH_LOC_4)
      else
        ifneq ($(strip $(ARCH_LOC_5)),)
          $(info Using $(ARCH_LOC_5))
          include $(ARCH_LOC_5)
        else
          $(error Could not find Makefile.arch!)
        endif
                endif
    endif
  endif
endif

DIR = .

BINDIR       = $(DIR)/bin
TMPDIR       = $(DIR)/tmp

ROOTLIBS     = $(shell $(ROOTSYS)/bin/root-config --libs)

OBJS         = $(TMPDIR)/Plots.o
OBJS         += $(TMPDIR)/SteerPlotter.o
OBJS         += $(TMPDIR)/SteerParser.o
OBJS         += $(TMPDIR)/BaseSteer.o
OBJS         += $(TMPDIR)/SteerPlotter_Dict.o
OBJS         += $(TMPDIR)/BaseSteer_Dict.o

OBJS         += $(TMPDIR)/SHist.o
OBJS         += $(TMPDIR)/SHist_Dict.o
OBJS         += $(TMPDIR)/SPlotter.o
OBJS         += $(TMPDIR)/FileParser.o

DICTFILE     = SteerPlotter_Dict.cxx
DICTFILE     += BaseSteer_Dict.cxx
DICTFILE     += SHist_Dict.cxx


all: setup Plots

setup:
	mkdir -p $(BINDIR)
	mkdir -p $(TMPDIR)

clean:
	@rm -f *Dict.cxx *Dict.h
	@rm -rf $(OBJS) $(BINDIR) $(TMPDIR)

distclean: clean

$(TMPDIR)/%.o: %.cxx
	$(CXX) -c $(CXXFLAGS) -I$(TMPDIR) \
	-o $@ $<

# Rule to create the dictionary
$(DICTFILE): $(patsubst %_Dict.cxx,%.h,$@)
	@echo "Generating dictionary $@ with file $(patsubst %_Dict.cxx,%.h,$@)"
	@$(shell root-config --exec-prefix)/bin/rootcint -f $@ -c -p -I/. $(patsubst %_Dict.cxx,%.h,$@)

Plots:  $(OBJS) $(DICTOBJ)
	$(CXX) $(CXXFLAGS) $(OBJS) $(DICTOBJ) $(ROOTLIBS) \
	-o $(BINDIR)/Plots
	@echo "Plotter for the B2G-16-015 created."
