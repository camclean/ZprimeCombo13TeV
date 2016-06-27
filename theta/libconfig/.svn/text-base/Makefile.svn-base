include ../Makefile.options

sources := $(wildcard *.c)
objects := $(patsubst %.c,.bin/%.o,$(sources)) .bin/libconfigcpp.o
deps := $(patsubst %.o,%.d,$(objects))

include ../Makefile.rules

#use fewer warning here:
WFLAGSCPP := -Wall

all: ../lib/libconfig.so

../lib/libconfig.so: $(objects)
	@$(LINK_SO)

clean-subdir:
	@rm -f ../lib/libconfig.so

