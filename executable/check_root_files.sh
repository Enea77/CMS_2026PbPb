#!/bin/bash

# Ensure a directory was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory_path>"
    exit 1
fi

DIR="$1"
OUTPUT_FILE="valid_files.txt"

# Ensure the directory exists
if [ ! -d "$DIR" ]; then
    echo "Error: Directory '$DIR' does not exist."
    exit 1
fi

# Clear or create the output file
> "$OUTPUT_FILE"

echo "Scanning directory: $DIR"
echo "Valid files will be saved to: $OUTPUT_FILE"
echo "---------------------------------------------------"

# Find all .root files and process them line by line
find "$DIR" -type f -name "*.root" | while read -r file; do
    
    # 1. Quick system check: skip 0-byte files immediately
    if [ ! -s "$file" ]; then
        echo "[EMPTY/0-BYTE]  $file"
        continue
    fi

    # 2. Deep ROOT check: Zombie, Recovered, or Zero Keys
    # We exit with status 1 if it's bad, and 0 if it's good.
    root -l -b -q -e "
        TFile *f = TFile::Open(\"$file\"); 
        if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered) || f->GetNkeys() == 0) { 
            gSystem->Exit(1); 
        }
    " > /dev/null 2>&1
    
    # Check the exit status of the ROOT command
    if [ $? -eq 0 ]; then
        # File is perfectly fine
        echo "$file" >> "$OUTPUT_FILE"
    else
        echo "[CORRUPT/NO KEYS] $file"
    fi
done

echo "---------------------------------------------------"
echo "Done! Found $(wc -l < "$OUTPUT_FILE") valid files."