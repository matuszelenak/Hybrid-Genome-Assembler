#!/usr/bin/env bash

if [ ! -f "$2" ]; then
  # Find alignments between pairs of reads
  minimap2 -x ava-"$3" -t 8 "$1" "$1" 2>/dev/null | gzip -1 > nanopore.paf.gz
  # Use overlaps to compute the assembled genome
  miniasm -f "$1" nanopore.paf.gz > miniasm.gfa 2> miniasm.log
  # Convert genome to fasta format
  perl -lane 'print ">$F[1]\n$F[2]" if $F[0] eq "S"' miniasm.gfa > "$2"
  ## Align reads to the assembled genome
  #minimap2 -x map-ont --secondary=no -t 8 "$2" "$1" | gzip -1 > miniasm.paf.gz
  ## Polish the genome by finding consensus of aligned reads at each position
  #racon -t 8 "$1" miniasm.paf.gz "$2" > "$2"

  rm nanopore.paf.gz miniasm.paf.gz  miniasm.gfa miniasm.log
fi

