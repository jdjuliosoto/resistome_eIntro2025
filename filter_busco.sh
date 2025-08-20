#!/bin/bash

BUSCO_DIR=$1
OUTFILE=$2

echo -e "bin_id\tC\tS\tD\tF\tM" > "$OUTFILE"

# Find all summary files within subfolders
find "$BUSCO_DIR" -type f -name "short_summary.specific.fungi.odb12.*.txt" | while read -r FILE; do
  BIN_ID=$(basename "$FILE" | sed 's/short_summary\.specific\.fungi\.odb12\.//; s/\.txt//')

  # Extract the summary line
  LINE=$(grep -Eo "C:[^ ]+" "$FILE")

  # Parsing values
  C=$(echo "$LINE" | sed -E 's/C:([0-9.]+)\[.*/\1/')
  S=$(echo "$LINE" | sed -E 's/.*S:([0-9.]+)%.*/\1/')
  D=$(echo "$LINE" | sed -E 's/.*D:([0-9.]+)%.*/\1/')
  F=$(grep -Eo "F:[0-9.]+%" "$FILE" | cut -d':' -f2 | tr -d '%')
  M=$(grep -Eo "M:[0-9.]+%" "$FILE" | cut -d':' -f2 | tr -d '%')

  echo -e "$BIN_ID\t$C\t$S\t$D\t$F\t$M" >> "$OUTFILE"
done