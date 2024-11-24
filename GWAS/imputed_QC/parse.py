#!/usr/bin/env python3

import gzip
import sys
from argparse import ArgumentParser

# Argument parsing
parser = ArgumentParser(description='Filter genotype file according to MAF and R2.')
parser.add_argument("-g", "--genotypes", type=str, help="Genotypes in VCF")
parser.add_argument("-m", "--maf", type=float, help="MAF threshold")
parser.add_argument("-r", "--r_square", type=float, help="R2 threshold")
parser.add_argument("-o", "--out", type=str, help="Output VCF")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Open input VCF
if args.genotypes.endswith('.gz'):
    fh = gzip.open(args.genotypes, 'rt')
else:
    fh = open(args.genotypes, 'r')

# Open output VCF
out_fh = open(args.out, 'w')

# Read line by line
for line in fh:
    if line.startswith('#'):
        # Print meta-information lines and header
        out_fh.write(line)
    else:
        # Obtain INFO field
        info = line.strip().split('\t')[7]
        # Filter MAF/R2. In my test VCF: MAF=0.1;R2=0.1;NS=3;DP=14
        infos = info.split(';')
        maf = float(infos[1].split("=")[1])
        r2 = float(infos[2].split("=")[1])
        if maf > args.maf and r2 > args.r_square:
            out_fh.write(line)

# Close files
fh.close()
out_fh.close()
