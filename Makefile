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

.PHONY: update_pages
update_pages: clean docs
	rm -rf stash
	mkdir -p stash/$(BUILD_ROOT)
	cp $(BUILD_ROOT)/cov stash/$(BUILD_ROOT)/
	cp $(BUILD_ROOT)/coverage.svg stash/$(BUILD_ROOT)/
	cp docs stash/
	git checkout gh-pages
	cp stash/$(BUILD_ROOT)/cov $(BUILD_ROOT)/
	cp stash/$(BUILD_ROOT)/coverage.svg $(BUILD_ROOT)/
	cp stash/docs ./

.PHONY: debug_makefile
debug_makefile:
	$(MAKE) -C $(BUILD_ROOT) $(MAKECMDGOALS)

