include ../Makefile.options

libsources := $(wildcard *.cpp)
libobjects := $(patsubst %.cpp,.bin/%.o,$(libsources))
deps := $(patsubst %.cpp,.bin/%.d,$(libsources))

libobjects += .bin/atomic.o

all: ../lib/libtheta.so

include ../Makefile.rules

# switch between implementations of atomic:
ifeq ($(OPTIONS_GENERIC_ARCH),yes)

.bin/atomic.o: atomic.cxx
	@$(COMPILE) $(FLAGSCPP) -o $@ $<

else

.bin/atomic.o: atomic.s
	@$(COMPILE) $(FLAGSC) -o $@ $<
.bin/atomic_asserts.o: atomic_asserts.cxx
	@$(COMPILE) $(FLAGSCPP) -o $@ $<

libobjects += .bin/atomic_asserts.o

endif



../lib/libtheta.so: $(libobjects)
	@$(LINK_SO) -lconfig


clean-subdir:
	@rm -f ../lib/libtheta.so

