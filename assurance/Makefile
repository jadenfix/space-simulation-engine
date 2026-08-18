CC ?= cc
CFLAGS ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror
CPPFLAGS ?= -Iinclude
LDFLAGS ?=
LDLIBS ?= -lm
BUILD := build
OUTPUT := output
LIBOBJ := $(BUILD)/core.o $(BUILD)/physics.o
TEST := $(BUILD)/test_assurance

.PHONY: all test clean asan ubsan
all: $(TEST)
$(BUILD) $(OUTPUT):
	mkdir -p $@
$(BUILD)/core.o: src/core.c include/spacewind/assurance.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
$(BUILD)/physics.o: src/physics.c include/spacewind/assurance.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
$(TEST): tests/test_assurance.c $(LIBOBJ) include/spacewind/assurance.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) tests/test_assurance.c $(LIBOBJ) $(LDLIBS) -o $@
test: $(TEST) | $(OUTPUT)
	./$(TEST)
	python3 -m json.tool $(OUTPUT)/assurance_receipt.json >/dev/null
	python3 -m json.tool $(OUTPUT)/cycle_receipt.json >/dev/null
asan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-std=c11 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror -fsanitize=address,undefined -fno-omit-frame-pointer' LDFLAGS='-fsanitize=address,undefined' test
ubsan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-std=c11 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror -fsanitize=undefined -fno-omit-frame-pointer' LDFLAGS='-fsanitize=undefined' test
clean:
	rm -rf $(BUILD) $(OUTPUT)/*.json $(OUTPUT)/*.csv
