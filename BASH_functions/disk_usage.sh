#!/usr/bin/env bash
# Usage: ./disk_usage.sh [path]
# Default path is current directory

TARGET="${1:-.}"

# Get absolute path and folder name
ABS_TARGET=$(realpath "$TARGET")
FOLDER_NAME=$(basename "$ABS_TARGET")

# Timestamp
DATE=$(date +"%Y-%m-%d_%H-%M-%S")

# Output file
OUTFILE="disk_usage_${FOLDER_NAME}_${DATE}.txt"

progress() {
    local count=0
    local spin='-\|/'
    local i=0
    while IFS= read -r line; do
        count=$((count+1))
        i=$(( (i+1) %4 ))
        last=$(echo "$line" | cut -f2-)
        printf "\r[%s] %d processed | Last: %s" \
            "${spin:$i:1}" "$count" "$last"
        echo "$line" >> "$1"
    done
    printf "\rDone. Processed %d entries. Output saved to %s\n" "$count" "$1"
}

# Header
{
    echo "### Disk usage report for: $ABS_TARGET"
    echo "### Generated: $(date)"
    echo ""
    echo "### Full recursive usage (files + folders) ###"
} > "$OUTFILE"

# Files + folders
du -ah "$TARGET" 2>/dev/null | sort -hr | progress "$OUTFILE"

# Folder-only ranking
{
    echo ""
    echo "### Folder size ranking (directories only) ###"
} >> "$OUTFILE"

du -h "$TARGET" 2>/dev/null | sort -hr | progress "$OUTFILE"

