#!/usr/bin/env bash

jellyfish bc -C -m "$2" -s 1G -t 8 -o temp.bc "$1"
jellyfish count -C -m "$2" -s 1G -t 8 --bc temp.bc -o mers.jf "$1"
jellyfish dump -c mers.jf -o mers.txt

rm temp.bc
rm mers.jf
