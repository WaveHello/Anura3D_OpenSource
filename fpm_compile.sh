#!/bin/bash

# This shell file makes it easier to run fpm and pipe the compile output to a file called fpm_compile_output.txt
# The purpose of fpm_compile_output.txt is to hold the possible warnings and errors to make it easier to search

FILE="@file"  # The file with the fpm compile instructions

# Run fpm with the specified file
fpm "$FILE" | tee fpm_compile_output.txt 

# Check if fpm succeeded
if [ $? -eq 0 ]; then
    echo "Compilation succeeded for $FILE"
else
    echo "Compilation failed for $FILE"
    exit 1
fi
