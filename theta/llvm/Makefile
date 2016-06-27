libsources := $(wildcard *.cpp)
libobjects := $(patsubst %.cpp,.bin/%.o,$(libsources))
deps := $(patsubst %.cpp,.bin/%.d,$(libsources))

include ../Makefile.options
include ../Makefile.rules

EXTRA_COMPILE += $(shell llvm-config-3.1 --cppflags)

all: ../lib/llvm-plugins.so

../lib/llvm-plugins.so: $(libobjects)
	@$(LINK_SO) $(shell llvm-config-3.1 --libfiles all) #nativecodegen engine scalaropts ipo)

test: .bin/test.o ../lib/llvm-plugins.so
	@$(LINK_EXE) ../lib/core-plugins.so

clean-subdir:
	@rm -f ../lib/llvm-plugins.so
