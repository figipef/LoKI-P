# Compiler
CXX := g++
CXXFLAGS := -std=c++14 -Wall -Wextra -Wpedantic -Wconversion -Wunused-parameter -Wunused-but-set-parameter -O2

# Directories
ROOTDIR := ..
BUILDDIR := $(ROOTDIR)/build
SRCDIR := $(ROOTDIR)/src
BINDIR := $(ROOTDIR)/bin
INCDIR := $(ROOTDIR)/inc
OBJDIR := $(BUILDDIR)/obj

# Files
SOURCES := $(wildcard $(SRCDIR)/*.cpp) #$(wildcard $(INCDIR)\chi2/*.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

$(info $(SRCDIR))

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
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

# Create build directory
$(OBJDIR):
	mkdir obj

# Create bin directory
$(BINDIR):
	@if not exist "$(BINDIR)" mkdir "$(BINDIR)"

# Clean
clean:
	rm -rf $(OBJDIR) $(TARGET)

# Phony targets
.PHONY: all clean
