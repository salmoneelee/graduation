# Compiler 
CC = clang

# Compiler Flags
CFLAGS = -g -Wall -Wextra -fsanitize=address -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lgsl -lgslcblas -lm -fsanitize=address

# Target Executable
TARGET = trial9

# Source and Header Files
SRC = trial9.c trial9_functions.c
HDR = trial9.h

# Files to Clean
CLEAN_FILES = position_data_9.txt text_data_9.txt flux_data_9.txt $(TARGET) 

# Build Target
$(TARGET): $(SRC) $(HDR)
	$(CC) -o $(TARGET) $(SRC) $(CFLAGS) $(LDFLAGS)

# Clean Up
clean:
	rm -f $(CLEAN_FILES)
