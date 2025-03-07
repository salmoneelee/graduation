# Compiler 
CC = clang

# Compiler Flags
CFLAGS = -g -Wall -Wextra -fsanitize=address -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lgsl -lgslcblas -lm -fsanitize=address

# Target Executable
TARGET = trial1 


# Source and Header Files
SRC = trial1.c trial1_functions.c
HDR = trial1.h

# Files to Clean
CLEAN_FILES = position_data.txt text_data.txt flux_data.txt $(TARGET) 

# Build Target
$(TARGET): $(SRC) $(HDR)
	$(CC) -o $(TARGET) $(SRC) $(CFLAGS) $(LDFLAGS)

# Clean Up
clean:
	rm -f $(CLEAN_FILES)
