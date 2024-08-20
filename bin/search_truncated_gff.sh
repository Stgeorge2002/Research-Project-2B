#!/bin/bash
set -e
echo "Starting search process" > search_log.txt
echo "Contents of truncated_genes file:" >> search_log.txt
cat $1 >> search_log.txt
echo "List of cleaned FASTA files:" >> search_log.txt
ls -l *.fa >> search_log.txt

while IFS= read -r genome_id || [[ -n "$genome_id" ]]; do
    # Extract the core genome ID pattern (e.g., "11657_2#28")
    core_id=$(echo "$genome_id" | grep -oP '^\d+_\d+#\d+')
    echo "Searching for genome_id: $core_id" >> search_log.txt
    
    # Use find in the current directory, following symlinks
    matching_file=$(find -L . -maxdepth 1 -type f -name "${core_id}*.fa" | head -n 1)
    
    if [ -n "$matching_file" ]; then
        new_name="truncated_$(basename $matching_file)"
        cp -L "$matching_file" "$new_name"
        echo "Copied $matching_file to $new_name" >> search_log.txt
    else
        echo "No matching file found for $core_id" >> search_log.txt
    fi
done < $1

echo "Finished search process" >> search_log.txt
echo "Files in current directory:" >> search_log.txt
ls -l >> search_log.txt