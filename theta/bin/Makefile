include ../Makefile.options

deps := .bin/merge.d .bin/theta.d

include ../Makefile.rules
-include Makefile.local

all: theta


merge: .bin/merge.o
	@$(LINK_EXE)

theta: .bin/theta.o
	@$(LINK_EXE)

clean-subdir:
	@rm -f merge theta

