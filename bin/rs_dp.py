#!/usr/bin/env python
import pysam
import argparse
from numpy import median

USAGE = """
rs_dp.py -i <input BAM file> -b <GBSBED file of restriction sites> -n <normfactors>
"""

DESCRIPTION = """
Parse reads in a BAM file and calculate depth of coverage at restriction sites
in a GBSBED file. Also calculate normalized depth of coverage according to the
normalization factors provided by the user.
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('-i', '--input', dest='i', help='input BAM file')
parser.add_argument('-s', '--sample', dest='sample', help='sample name')
parser.add_argument('-b', '--bed', dest='bed', help='GBSBED file of restriction sites')
parser.add_argument('-n', '--normfactors', dest='nf', help='DP normalization factors')
parser.add_argument('-w', '--window', dest='window', type=int, default=102, help='window for searching for reads around restriction sites (default=102)')
args = parser.parse_args()

# Parse the normfactors file and store in dict.
normfactors = {}
nf = open(args.nf, 'r')
header = nf.readline()
header = header.strip('#').strip()
header = header.split()
sample_index = header.index(args.sample)
for line in nf:
    line = line.strip()
    fields = line.split()
    insert = int(fields[0])
    normfactors[insert] = float(fields[sample_index])

header = ['chrom', 'start', 'end', 'fwd_insert', 'rev_insert', 'enzyme', 'strand',
          'fwd_ins_med', 'fwd_counts', 'fwd_normcounts',
          'rev_ins_med', 'rev_counts', 'rev_normcounts']
print '#' + '\t'.join(header)
window = args.window
bed = open(args.bed, 'r')
sam = pysam.Samfile(args.i, 'rb')
for line in bed:
    line = line.strip()
    # Ignore header.
    if line[0] == "#":
        continue
    (chrom, start, end, fwd_ins, rev_ins, enzyme,
     strand, site_id, freq, cut_allele, fwd_lig, rev_lig) = line.split()
    start = int(start)
    end = int(end)
    site_tag = "%s;%s;%s" % (enzyme, strand, start)
    # Read counts, keyed by read.is_reverse.
    counts = {0:0, 1:0}
    # Insert lists, keyed by read.is_reverse.
    inserts = {0:[], 1:[]}
    for read in sam.fetch(chrom, max(start - window, 0), end + window):
        tags = dict(read.tags)
        if 'Z1' in tags:
            print tags
        try:
            read_tag = tags['Z2']
        except:
            read_tag = None
        try:
            insert = tags['Z0']
        except:
            insert = None
        if read_tag == site_tag:
            counts[read.is_reverse] += 1
            inserts[read.is_reverse].append(insert)
    if inserts[0]:
        fwd_insert_med = median([i for i in inserts[0] if i])
    else:
        fwd_insert_med = '.'
    if inserts[1]:
        rev_insert_med = median([i for i in inserts[1] if i])
    else:
        rev_insert_med = '.'

    try:
        fwd_normcounts = counts[0] / normfactors[int(fwd_insert_med)]
        fwd_normcounts = '{0:.3f}'.format(fwd_normcounts)
    except:
        fwd_normcounts = '.'
    try:
        rev_normcounts = counts[1] / normfactors[int(rev_insert_med)]
        rev_normcounts = '{0:.3f}'.format(rev_normcounts)
    except:
        rev_normcounts = '.'
    output = (chrom, start, end, fwd_ins, rev_ins, enzyme, strand, 
              fwd_insert_med, counts[0], fwd_normcounts,
              rev_insert_med, counts[1], rev_normcounts)
    print '\t'.join([str(i) for i in output])
