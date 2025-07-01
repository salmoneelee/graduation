# Compiler 
CC = clang

# Compiler Flags
CFLAGS = -g -Wall -Wextra -fsanitize=address -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lgsl -lgslcblas -lm -fsanitize=address

# Target Executable
TARGET = trial

# Source and Header Files
SRC = trial.c trial_functions.c
HDR = trial.h

# Files to Clean
CLEAN_FILES = position_data.txt text_data.txt flux_data.txt k_data.txt reweight_bank.txt $(TARGET) 

# Build Target
$(TARGET): $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Clean Up
clean:
	rm -f $(CLEAN_FILES)
