# Compiler 
CC = clang

# Compiler Flags
CFLAGS = -g -Wall -Wextra -fsanitize=address -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lgsl -lgslcblas -lm -fsanitize=address

# Target Executable
TARGET = trial4

# Source and Header Files
SRC = trial4.c trial4_functions.c
HDR = trial4.h

# Files to Clean
CLEAN_FILES = position_data_4.txt text_data_4.txt flux_data_4.txt $(TARGET) 

# Build Target
$(TARGET): $(SRC) $(HDR)
	$(CC) -o $(TARGET) $(SRC) $(CFLAGS) $(LDFLAGS)

# Clean Up
clean:
	rm -f $(CLEAN_FILES)
