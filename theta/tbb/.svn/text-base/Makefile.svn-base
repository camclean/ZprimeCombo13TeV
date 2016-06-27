libsources := $(wildcard *.cpp)
libobjects := $(patsubst %.cpp,.bin/%.o,$(libsources))
deps := $(patsubst %.cpp,.bin/%.d,$(libsources))

include ../Makefile.options
include ../Makefile.rules

EXTRA_COMPILE := -std=c++0x

all: ../lib/tbb-plugins.so

../lib/tbb-plugins.so: $(libobjects)
	@$(LINK_SO) -ltbb

clean-subdir:
	@rm -f ../lib/tbb-plugins.so
