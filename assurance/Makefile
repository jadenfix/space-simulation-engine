CC ?= cc
CFLAGS ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror
CPPFLAGS ?= -Iinclude
LDFLAGS ?=
LDLIBS ?= -lm

BUILD := build
OUTPUT := output
LIBSRC := \
	src/core.c \
	src/physics.c \
	src/invariants.c \
	src/experiment.c \
	src/evidence.c \
	src/field_ledger.c \
	src/charge_deposition.c
LIBOBJ := $(patsubst src/%.c,$(BUILD)/%.o,$(LIBSRC))
TEST_NAMES := assurance invariants experiment evidence field_ledger charge_deposition
TEST_BINS := $(addprefix $(BUILD)/test_,$(TEST_NAMES))

.PHONY: all test clean asan ubsan
all: $(TEST_BINS)

$(BUILD) $(OUTPUT):
	mkdir -p $@

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_%: tests/test_%.c $(LIBOBJ) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< $(LIBOBJ) $(LDLIBS) -o $@

test: $(TEST_BINS) | $(OUTPUT)
	@set -e; for test_bin in $(TEST_BINS); do ./$$test_bin; done
	@set -e; for receipt in $(OUTPUT)/*.json; do python3 -m json.tool "$$receipt" >/dev/null; done

asan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-std=c11 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror -fsanitize=address,undefined -fno-omit-frame-pointer' LDFLAGS='-fsanitize=address,undefined' test

ubsan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-std=c11 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror -fsanitize=undefined -fno-omit-frame-pointer' LDFLAGS='-fsanitize=undefined' test

clean:
	rm -rf $(BUILD) $(OUTPUT)/*.json $(OUTPUT)/*.csv
