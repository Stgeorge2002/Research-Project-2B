#!/bin/bash
set -e
set -o pipefail

# Use the provided search terms and retmax
SEARCH_TERM="$1"
ACETYLTRANSFERASE_SEARCH_TERM="$2"
RETMAX="$3"
ENCODED_TERM=$(echo "$ACETYLTRANSFERASE_SEARCH_TERM" | sed 's/ /%20/g; s/\[/%5B/g; s/\]/%5D/g')

# Perform the initial esearch to get WebEnv and QueryKey
esearch_result=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=${ENCODED_TERM}&usehistory=y&retmax=${RETMAX}")
WEBENV=$(echo "$esearch_result" | grep -oPm1 "(?<=<WebEnv>)[^<]+")
QUERYKEY=$(echo "$esearch_result" | grep -oPm1 "(?<=<QueryKey>)[^<]+")

if [ -z "$WEBENV" ] || [ -z "$QUERYKEY" ]; then
    echo "Error: Failed to retrieve WebEnv or QueryKey"
    exit 1
fi

# Fetch the sequences
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&WebEnv=${WEBENV}&query_key=${QUERYKEY}&rettype=fasta_cds_na&retmode=text&retmax=${RETMAX}" > raw_sequences.fasta

# Filter the sequences
awk -v search_term="$SEARCH_TERM" '/^>/ {if (p) {print s; s=""}; p=index(tolower($0), tolower(search_term)); s=$0; next} p {s=s"\n"$0} END {if (p) print s}' raw_sequences.fasta > filtered_sequences.fasta

# Create a tab-delimited file with gene headers and sequences
awk '/^>/ {if (seq) print header "\t" seq; header = substr($0,2); seq=""; next} {seq = seq $0} END {print header "\t" seq}' filtered_sequences.fasta > gene_info.tsv

# Verify the results
total_seqs=$(grep -c "^>" filtered_sequences.fasta)
matching_seqs=$(grep -ic "$SEARCH_TERM" filtered_sequences.fasta)

echo "Total sequences: $total_seqs"
echo "Sequences with '$SEARCH_TERM' in protein field: $matching_seqs"

if [ "$total_seqs" -eq "$matching_seqs" ]; then
    echo "All retrieved sequences contain '$SEARCH_TERM' in the protein field."
else
    echo "Error: Not all sequences contain '$SEARCH_TERM' in the protein field."
    exit 1
fi

# Create BLAST nucleotide database
makeblastdb -in filtered_sequences.fasta \
            -dbtype nucl \
            -out spneumo_custom_db \
            -title "S. pneumoniae custom nucleotide DB"

# Check if the database was created successfully
if [ ! -f spneumo_custom_db.nhr ]; then
    echo "Error: Failed to create BLAST database"
    exit 1
fi

echo "Custom S. pneumoniae nucleotide database created successfully"