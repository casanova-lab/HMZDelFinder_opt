#! /usr/bin/env bash

# Convert a list of files to YAML format:
# ID: path/to/file

while read path; do
    id=$(basename "$path" | sed 's/[-_\.].*//g')
    echo "${id}: $path"
done < "$1"
