include ./Makefile.options
-include ./Makefile.options.local

#order matters
DIRS = libconfig src plugins root bin test

ifeq ($(OPTIONS_BUILD_LLVM),yes)
   DIRS += llvm
endif

ifeq ($(OPTIONS_BUILD_TBB),yes)
   DIRS += tbb
endif

all:
	@for d in $(DIRS); do ( $(MAKE) -C $$d ) || break; done

clean:
	@for d in $(DIRS); do ( $(MAKE) -C $$d clean ) || true; done
	@cd test/test-stat; $(MAKE) clean

#ok, this is not very nice, as it hardcodes the documentation copy target, but as long as I am the only
# developer, it should work well:
doc:
	@doxygen
	@( warn=`wc -l doxygen-warn.txt | cut -f 1 -d" "`; if [ $$warn -gt 0 ]; then echo There have been about $$warn warnings from doxygen, see doxygen-warn.txt; fi )

#	@cp doc/tabs.css doc/html
#	@if [ "`hostname`" = "ekplx22" ]; then \
#	   rsync -a --del doc/* /usr/users/ott/public_html/theta/testing; \
#	fi
   

.PHONY: clean all doc

