#!/usr/bin/python

import sys
import os
import optparse
import GFFutils
usage = """

    Clean a GFF file by doing the following:

        - remove FASTA sequence at the end
        - remove comments
        - remove items that have start > stop
        - optionally remove features with negative start position
        - optionally remove items from a list of feature types to exclude
        - optionally add "chr" to the beginning of chromosomes

    Default output is to stdout unless the -o option is specified.

    Typical usage:

        %s --addchr --exclude orthologous_to,pcr_product,BAC_cloned_genomic_insert -o cleaned.gff input.gff

""" % sys.argv[0]
op = optparse.OptionParser(usage=usage)
op.add_option('--addchr', action='store_true', help='Prefix each chromosome with "chr" '
                                                    'in the output file.')
op.add_option('-o',help='Optional output file; if unspecified output will go to stdout.')
op.add_option('--exclude', help='Comma-separated list of feature types '
                                'to exclude from the output.')
op.add_option('--no-negative', dest='no_neg', action='store_true', help='Remove features that have negative start or stop positions')

options,args = op.parse_args()
if len(args) != 1:
    op.print_help()
    sys.exit(1)

gfffn = args[0]

if options.o is None:
    out = sys.stdout
else:
    out = open(options.o,'w')


# list of featuretypes that you want to remove.
culled_features = options.exclude.split(',')

f = GFFutils.GFFFile(gfffn)

for feature in f:
    if feature.start is None:
        continue
    if feature.stop is None:
        continue
    if feature.start > feature.stop:
        continue
    if feature.featuretype in culled_features:
        continue
    if options.addchr:
        feature.chr = 'chr'+feature.chr
    if options.no_neg:
        if (feature.start < 0) or (feature.stop < 0):
            continue
    out.write(feature.tostring())
    out.flush()
 
if options.o is not None:
    out.close()
