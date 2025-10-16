# Compiler
CXX := mpicxx
CXXFLAGS := -std=c++14 -Wall -Wextra -Wpedantic -Wconversion -Wunused-parameter -Wunused-but-set-parameter -O2 \
            -I../inc

# Directories
ROOTDIR := ..
BUILDDIR := $(ROOTDIR)/build
SRCDIR := $(ROOTDIR)/src
BINDIR := $(ROOTDIR)/bin
OBJDIR := $(BUILDDIR)/obj

# Files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Program name
PROGRAM_NAME := test
ifdef NAME
    PROGRAM_NAME := $(NAME)
endif

# Target binary
TARGET := $(BINDIR)/$(PROGRAM_NAME)

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CXX) $(OBJECTS) -o $@

# Compilation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# Clean
clean:
	rm -rf $(OBJDIR) $(TARGET)

# Phony targets
.PHONY: all clean
