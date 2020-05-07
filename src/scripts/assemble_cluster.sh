#!/usr/bin/env bash

if [ ! -f "$2" ]; then
  # Find alignments between pairs of reads
  minimap2 -x ava-pb -t 8 "$1" "$1" 2>/dev/null | gzip -1 > mapping.paf.gz
  # Use overlaps to compute the assembled genome
  miniasm -f "$1" mapping.paf.gz > miniasm.gfa 2> miniasm.log
  # Convert genome to fasta format
  perl -lane 'print ">$F[1]\n$F[2]" if $F[0] eq "S"' miniasm.gfa > unpolished.fa
  # Align reads to the assembled genome
  minimap2 -x map-pb --secondary=no -t 8 unpolished.fa "$1" | gzip -1 > mapping.paf.gz
  # Polish the genome by finding consensus of aligned reads at each position
  racon -t 8 "$1" mapping.paf.gz unpolished.fa > "$2"

  rm mapping.paf.gz miniasm.gfa miniasm.log unpolished.fa
fi
