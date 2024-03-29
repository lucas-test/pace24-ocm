#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "Error: No directory specified."
    echo "Usage: ./all.sh <directory>"
    exit 1
fi

> results.out

# dir=$1
# for file in "$dir"/*.gr; do
#     ./bin/main "$file" >> results.out
# done

timeout_duration=30
dir=$1
for file in `ls "$dir" | sort -g`; do
    SECONDS=0
    timeout $timeout_duration ./bin/main "$dir""$file" >> results.out
    if [ $? -eq 0 ]; then
        echo "Execution of ./bin/main for $file took $SECONDS seconds." >> results.out
    elif [ $? -eq 124 ]; then # 124 is the timeout exit status
        echo "Execution of ./bin/main for $file exceeded $timeout_duration seconds." >> results.out
    fi
done

