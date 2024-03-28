#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "Error: No directory specified."
    echo "Usage: ./all.sh <directory>"
    exit 1
fi

> results.out

dir=$1
for file in "$dir"/*.gr; do
    ./bin/main "$file" >> results.out
done