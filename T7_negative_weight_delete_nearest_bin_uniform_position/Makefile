# Compiler 
CC = clang

# Compiler Flags
CFLAGS = -g -Wall -Wextra -fsanitize=address -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lgsl -lgslcblas -lm -fsanitize=address

# Target Executable
TARGET = trial7

# Source and Header Files
SRC = trial7.c trial7_functions.c
HDR = trial7.h

# Files to Clean
CLEAN_FILES = position_data_7.txt text_data_7.txt flux_data_7.txt $(TARGET) 

# Build Target
$(TARGET): $(SRC) $(HDR)
	$(CC) -o $(TARGET) $(SRC) $(CFLAGS) $(LDFLAGS)

# Clean Up
clean:
	rm -f $(CLEAN_FILES)
