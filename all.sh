#!/bin/bash

# Initialize an empty file to store the results
> results.out

# Loop through each .gr file in the graphs/exact directory
for file in graphs/exact/*.gr; do
    # Execute the command for each file and append the output to results.out
    ./bin/main "$file" >> results.out
done
