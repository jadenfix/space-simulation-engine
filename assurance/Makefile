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
	src/charge_deposition.c \
	src/spectral_poisson.c \
	src/maxwell1d_reference.c \
	src/covariant_em.c \
	src/em_pic1d.c
LIBOBJ := $(patsubst src/%.c,$(BUILD)/%.o,$(LIBSRC))
TEST_NAMES := assurance invariants experiment evidence field_ledger charge_deposition spectral_poisson maxwell1d_reference covariant_em maxwell_covariant_integration em_pic1d em_pic1d_covariant_integration
UNIT_TEST_BINS := $(addprefix $(BUILD)/test_,$(TEST_NAMES))
PIC_TRANSPORT_TEST_BIN := $(BUILD)/test_pic_transport_integration
SPECTRAL_INTEGRATION_TEST_BIN := $(BUILD)/test_spectral_poisson_integration
FDTD_SPECTRAL_TEST_BIN := $(BUILD)/test_fdtd_spectral_integration
INTEGRATION_TEST_BINS := $(PIC_TRANSPORT_TEST_BIN) $(SPECTRAL_INTEGRATION_TEST_BIN) $(FDTD_SPECTRAL_TEST_BIN)
TEST_BINS := $(UNIT_TEST_BINS) $(INTEGRATION_TEST_BINS)

.PHONY: all test clean asan ubsan
all: $(TEST_BINS)

$(BUILD) $(OUTPUT):
	mkdir -p $@

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_%: tests/test_%.c $(LIBOBJ) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< $(LIBOBJ) $(LDLIBS) -o $@

$(PIC_TRANSPORT_TEST_BIN): tests/test_pic_transport_integration.c $(LIBOBJ) ../src/plasma.c ../src/math3.c | $(BUILD)
	$(CC) $(CPPFLAGS) -I../include $(CFLAGS) $(LDFLAGS) \
		tests/test_pic_transport_integration.c ../src/plasma.c ../src/math3.c \
		$(LIBOBJ) $(LDLIBS) -o $@

$(SPECTRAL_INTEGRATION_TEST_BIN): tests/test_spectral_poisson_integration.c $(LIBOBJ) ../src/plasma.c ../src/math3.c | $(BUILD)
	$(CC) $(CPPFLAGS) -I../include $(CFLAGS) $(LDFLAGS) \
		tests/test_spectral_poisson_integration.c ../src/plasma.c ../src/math3.c \
		$(LIBOBJ) $(LDLIBS) -o $@

$(FDTD_SPECTRAL_TEST_BIN): tests/test_fdtd_spectral_integration.c $(LIBOBJ) ../src/plasma.c ../src/math3.c | $(BUILD)
	$(CC) $(CPPFLAGS) -I../include $(CFLAGS) $(LDFLAGS) \
		tests/test_fdtd_spectral_integration.c ../src/plasma.c ../src/math3.c \
		$(LIBOBJ) $(LDLIBS) -o $@

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
