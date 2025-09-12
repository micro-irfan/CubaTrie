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

# ---- toolchain ----
CC      ?= gcc
CSTD     = -std=c11
WARN     = -Wall -Wextra -Werror=implicit-function-declaration
# (Optional) expose POSIX funcs like strndup/strnlen; uncomment if you use them:
# CPPFLAGS += -D_POSIX_C_SOURCE=200809L
CPPFLAGS += -MMD -MP              # auto deps

# ---- sanitizers (used in debug) ----
SAN      = -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all

# ---- modes ----
MODE ?= release
ifeq ($(MODE),debug)
  CFLAGS  += -g3 -O0 $(SAN) -D_POSIX_C_SOURCE=200809L
else ifeq ($(MODE),release)
  CFLAGS  += -O3 -g -D_POSIX_C_SOURCE=200809L -march=native -DNDEBUG
endif

CFLAGS  += $(CSTD) $(WARN)
LDFLAGS +=
LDLIBS  += -lz                     # zlib

# ---- rules ----
.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# compile C -> o (use same flags so ASan/UBSan link in too)
%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(DEPS) 

-include $(DEPS)
