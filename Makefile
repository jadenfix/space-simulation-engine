CC ?= cc
AR ?= ar
BUILD_DIR ?= build
CPPFLAGS ?= -Iinclude
CFLAGS ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes
LDFLAGS ?=
LDLIBS ?= -lm

LIB_SOURCES := \
	src/math3.c \
	src/rng.c \
	src/gravity_grid.c \
	src/nbody.c \
	src/metric.c \
	src/fields.c \
	src/spacecraft.c \
	src/integrator.c \
	src/plasma.c \
	src/mhd1d.c \
	src/simulation.c
LIB_OBJECTS := $(LIB_SOURCES:%.c=$(BUILD_DIR)/%.o)
CLI_OBJECT := $(BUILD_DIR)/src/cli.o
TEST_OBJECT := $(BUILD_DIR)/tests/test_main.o
LIBRARY := $(BUILD_DIR)/libspacewind.a
BINARY := $(BUILD_DIR)/spacewind
TEST_BINARY := $(BUILD_DIR)/spacewind_tests

.PHONY: all clean test asan ubsan examples validate validate-all package format-check

all: $(BINARY) $(TEST_BINARY)

$(BINARY): $(CLI_OBJECT) $(LIBRARY)
	$(CC) $(LDFLAGS) -o $@ $(CLI_OBJECT) $(LIBRARY) $(LDLIBS)

$(TEST_BINARY): $(TEST_OBJECT) $(LIBRARY)
	$(CC) $(LDFLAGS) -o $@ $(TEST_OBJECT) $(LIBRARY) $(LDLIBS)

$(LIBRARY): $(LIB_OBJECTS)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $(LIB_OBJECTS)

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

test: $(TEST_BINARY)
	$(TEST_BINARY)

asan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-std=c11 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -fsanitize=address,undefined -fno-omit-frame-pointer' LDFLAGS='-fsanitize=address,undefined' test

ubsan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-std=c11 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -fsanitize=undefined -fno-omit-frame-pointer' LDFLAGS='-fsanitize=undefined' test

examples: $(BINARY)
	@mkdir -p output
	$(BINARY) wind output/wind.csv
	$(BINARY) boris output/boris.csv
	$(BINARY) fdtd output/fdtd.csv
	$(BINARY) pic output/pic.csv
	$(BINARY) mhd output/mhd.csv
	$(BINARY) geodesic output/geodesic.csv
	$(BINARY) kerr output/kerr.csv
	$(BINARY) nbody output/nbody.csv
	$(BINARY) poisson output/gravity_grid.csv
	$(BINARY) simulate configs/baseline_1au.cfg output/baseline.csv output/baseline.receipt.json
	$(BINARY) simulate configs/hybrid_glider.cfg output/hybrid_glider.csv output/hybrid_glider.receipt.json

validate:
	./scripts/validate.sh

validate-all:
	./scripts/validate_all.sh

package:
	./scripts/package.sh

clean:
	rm -rf $(BUILD_DIR) output/*.csv output/*.json
