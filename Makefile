CC := gcc
LD := gcc
AR := ar
BUILD_ROOT := build/debug

SOURCES := $(shell find src -name '*.c')
HEADERS := $(shell find include -name '*.h')

.DEFAULT_GOAL := all

.PHONY: all
all:
	$(MAKE) -C $(BUILD_ROOT)

.PHONY: test
test:
	$(MAKE) -C $(BUILD_ROOT) $(MAKECMDGOALS)

.PHONY: coverage
coverage:
	$(MAKE) -C $(BUILD_ROOT) $(MAKECMDGOALS)

.PHONY: clean
clean:
	@for a in $$(ls build); do\
		if [ -d build/$$a ]; then \
			$(MAKE) -C build/$$a $(MAKECMDGOALS); \
		fi; \
	done;
	-rm -rf docs

docs: doxygen.conf $(SOURCES) $(HEADERS)
	doxygen $<
	mv html $@

.PHONY: debug_makefile
debug_makefile:
	$(MAKE) -C $(BUILD_ROOT) $(MAKECMDGOALS)

