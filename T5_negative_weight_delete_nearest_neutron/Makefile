# Compiler 
CC = clang

# Compiler Flags
CFLAGS = -g -Wall -Wextra -fsanitize=address -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lgsl -lgslcblas -lm -fsanitize=address

# Target Executable
TARGET = trial5

# Source and Header Files
SRC = trial5.c trial5_functions.c
HDR = trial5.h

# Files to Clean
CLEAN_FILES = position_data_5.txt text_data_5.txt flux_data_5.txt $(TARGET) 

# Build Target
$(TARGET): $(SRC) $(HDR)
	$(CC) -o $(TARGET) $(SRC) $(CFLAGS) $(LDFLAGS)

# Clean Up
clean:
	rm -f $(CLEAN_FILES)
