#!/usr/bin/env bash

minimap2 -x ava-pb -t 8 "$1" "$1" 2>/dev/null | gzip -1 > gzipfile

miniasm -f "$1" gzipfile > miniasm.gfa 2>/dev/null

perl -lane 'print ">$F[1]\n$F[2]" if $F[0] eq "S"' miniasm.gfa > "$2"

rm miniasm.gfa gzipfile