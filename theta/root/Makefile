include ../Makefile.options

ifeq ($(shell which root-config || echo -n A),A)
all:
	@echo -e \\nNOTE: no root-config found. Not building root-dependent stuff.\\n

include ../Makefile.rules

else

sources := $(wildcard *.cpp)
objects := $(patsubst %.cpp,.bin/%.o,$(sources))
deps := $(patsubst %.cpp,.bin/%.d,$(sources)) .bin/histos.d

include ../Makefile.rules

# just use -Wall here, do not warn about old-stryle casts, etc. Otherwise, we will be swamped by warnings
# from root.
WFLAGSCPP := -Wall
EXTRA_COMPILE := $(shell root-config --cflags)
EXTRA_LINK :=  -lMinuit2 -lMinuit $(shell root-config --libs) $(shell root-config --ldflags) 

# put the '-lNet' AFTER -lTree (this is required for the linker on oneiric which is sensitive to linking order):
EXTRA_LINK := $(patsubst -lNet,,$(EXTRA_LINK))
EXTRA_LINK := $(patsubst -lTree,-lTree -lNet,$(EXTRA_LINK))

all: ../lib/root.so

ifeq ($(OPTIONS_BUILD_TESTS),yes)
all: create_testhistos create_testtree
endif

../bin/histos: .bin/histos.o
	@$(LINK_EXE) $(EXTRA_LINK)

create_testhistos: .bin/create_testhistos.o
	@$(LINK_EXE_NOLIBS) $(EXTRA_LINK)

create_testhistos2: .bin/create_testhistos2.o
	@$(LINK_EXE_NOLIBS) $(EXTRA_LINK)

create_histos: .bin/create_histos.o
	@$(LINK_EXE_NOLIBS) $(EXTRA_LINK)

create_testtree: .bin/create_testtree.o
	@$(LINK_EXE_NOLIBS) $(EXTRA_LINK)

../lib/root.so: $(objects)
	@$(LINK_SO) -ltheta $(EXTRA_LINK)

endif

#define clean target in any case:
clean-subdir:
	@rm -f histos create_testhistos ../lib/root.so
