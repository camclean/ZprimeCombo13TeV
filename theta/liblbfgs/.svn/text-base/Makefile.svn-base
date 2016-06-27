objects := .bin/lbfgs.o .bin/lbfgs_minimizer.o
deps := $(patsubst %.o,%.d,$(objects))

include ../Makefile.rules

#all: 

all: ../lib/liblbfgs.so

../lib/liblbfgs.so: $(objects)
	@$(LINK_SO)

clean-subdir:
	@rm -f ../lib/liblbfgs.so
