#!/bin/bash

NODEFILE="${1:-PBS_NODEFILE}"

: > "$NODEFILE"

for nid in $(seq 1000 1006); do
    node=$(printf "nid%06d" "$nid")

    for ((slot=0; slot<192; slot++)); do
        echo "$node" >> "$NODEFILE"
    done
done

echo "Created $NODEFILE"
echo "Nodes: 7"
echo "PPN:   192"
echo "Slots: $(wc -l < "$NODEFILE")"