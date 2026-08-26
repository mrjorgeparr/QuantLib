# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -fopenmp
INCLUDES = -Iinclude
TARGET = libMyQuantLib.a
SRCDIR = src
INCDIR = include
OBJDIR = obj

# Source files
SOURCES = $(shell find $(SRCDIR) -name '*.cpp')
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Create directories if they don't exist
$(shell mkdir -p $(OBJDIR))

# Build target
$(TARGET): $(OBJECTS)
	ar rcs $@ $^

# Build object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean build files
clean:
	rm -rf $(OBJDIR) $(TARGET)

# Install dependencies (optional)
deps:
	@echo "No external dependencies required for this project"

# Debug build
debug: CXXFLAGS += -g -DDEBUG
debug: $(TARGET)

# Run the program
run: $(TARGET)
	@echo "Built $(TARGET). Link it from an executable or tests."

# Phony targets
.PHONY: clean deps debug run
