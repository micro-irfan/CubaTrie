# Makefile
# Usage:
#   make                # debug + sanitizers (default)
#   make MODE=release   # optimized build
#   make clean

# ---- project ----
TARGET := cuba_trie
SRCS   := kalloc.c utils.c trie.c kmer.c readseq.c main.c 
OBJS   := $(SRCS:.c=.o)
DEPS   := $(OBJS:.o=.d)
TEST_TARGET := cuba_trie_tests
TEST_SRCS   := tests/test_main.c trie.c kmer.c utils.c kalloc.c
GOLDEN_TEST_TARGET := cuba_trie_golden_tests
GOLDEN_TEST_SRCS   := tests/golden_test.c readseq.c trie.c kmer.c utils.c kalloc.c
TEST_DEPS := $(TEST_TARGET).d $(GOLDEN_TEST_TARGET).d

# ---- toolchain ----
CC      ?= gcc
CSTD     = -std=c11
WARN     = -Wall -Wextra -Werror=implicit-function-declaration
# (Optional) expose POSIX funcs like strndup/strnlen; uncomment if you use them:
# CPPFLAGS += -D_POSIX_C_SOURCE=200809L
CPPFLAGS += -MMD -MP              # auto deps

# ---- sanitizers (used in debug) ----
SAN      = -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all
THREAD_FLAGS = -pthread

# ---- modes ----
MODE ?= release
ifeq ($(MODE),debug)
  CFLAGS  += -g3 -O0 $(SAN) -D_POSIX_C_SOURCE=200809L
else ifeq ($(MODE),release)
  CFLAGS  += -O3 -g -D_POSIX_C_SOURCE=200809L -march=native -DNDEBUG
endif

CFLAGS  += $(CSTD) $(WARN) $(THREAD_FLAGS)
LDFLAGS += $(THREAD_FLAGS)
LDLIBS  += -lz                     # zlib

# ---- rules ----
.PHONY: all clean test golden-test
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# compile C -> o (use same flags so ASan/UBSan link in too)
%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(DEPS) $(TEST_TARGET) $(GOLDEN_TEST_TARGET) $(TEST_DEPS)

test: $(TEST_TARGET)
	@if [ -x "./$(TEST_TARGET)" ]; then \
		"./$(TEST_TARGET)"; \
	elif [ -x "./$(TEST_TARGET).exe" ]; then \
		"./$(TEST_TARGET).exe"; \
	else \
		echo "Test binary not found: $(TEST_TARGET) (or $(TEST_TARGET).exe)"; \
		exit 127; \
	fi

$(TEST_TARGET): $(TEST_SRCS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $(TEST_SRCS) $(LDLIBS)

golden-test: $(GOLDEN_TEST_TARGET)
	@if [ -x "./$(GOLDEN_TEST_TARGET)" ]; then \
		"./$(GOLDEN_TEST_TARGET)"; \
	elif [ -x "./$(GOLDEN_TEST_TARGET).exe" ]; then \
		"./$(GOLDEN_TEST_TARGET).exe"; \
	else \
		echo "Golden test binary not found: $(GOLDEN_TEST_TARGET) (or $(GOLDEN_TEST_TARGET).exe)"; \
		exit 127; \
	fi

$(GOLDEN_TEST_TARGET): $(GOLDEN_TEST_SRCS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $(GOLDEN_TEST_SRCS) $(LDLIBS)

-include $(DEPS) $(TEST_DEPS)
