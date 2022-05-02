#!/bin/env python3.10
import sys
import msprime
import tskit


#  ts = tskit.load(snakemake.input[0])
#  fh = open(str(snakemake.output[0]),"w")
#  tskit.write_ms(tree_sequence, snakemake.output[0])

ts = tskit.load(sys.argv[1])
fh = open(sys.argv[2],"w")
tskit.write_ms(ts, fh)


