#!/bin/bash

# Find crabs executable
CRABS_PATH=$(which crabs)
if [ -z "$CRABS_PATH" ]; then
    echo "Error: Could not locate crabs executable"
    exit 1
fi

# Create backup with timestamp
cp "$CRABS_PATH" "${CRABS_PATH}.backup.$(date +%Y%m%d_%H%M%S)"

# Apply the three changes:
# 1. Add whole_percent variable after the print statement
# 2. Replace first float(COV) comparison
# 3. Replace second float(COV) comparison
sed -i.bak '
    /print(f'\''filtering alignments based on parameter settings'\'')/ a\    whole_percent = float(COV)*100
    s/tcov >= float(COV) and/tcov >= whole_percent and/
    s/tcov >= float(COV)"/tcov >= whole_percent"/
' "$CRABS_PATH"

# Restore executable permissions
chmod +x "$CRABS_PATH"

echo "Patch applied successfully to $CRABS_PATH"