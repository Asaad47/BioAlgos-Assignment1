#!/bin/bash

srr_ids=(
  SRR11412973 SRR11412976 SRR11412979 SRR11412980 SRR11412984
  SRR21907296 SRR21907303 SRR21907307 SRR21907332 SRR21907330
)

for srr in "${srr_ids[@]}"; do
  echo "Fetching download links for $srr..."
  urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2 | tr ';' '\n')

  for url in $urls; do
    filename=$(basename "$url")
    echo "Downloading $filename ..."
    wget -c "ftp://${url}"
  done
done
