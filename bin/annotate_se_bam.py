#!/usr/bin/env python
import pysam
from sys import stdout
import argparse
from collections import namedtuple

# Parse command line arguments from user.
USAGE = """
annotate_se_bam.py -i <input bam file>
                   -o <output bam file>
                   -b <tabix-indexed restriction site bed file> (from make_gbsbed.py)
"""

DESCRIPTION = """
Query the mapping positions of reads in a BAM file against a tabix-indexed
GBSBED file of restriction sites (from make_gbsbed.py). Add the following
tags to each read (all positions are 0-indexed).

Z0: insert/fragment size (corrected for trimming)
Z1: ligation site (mapping position of the end base in the SAM ``Template``)
Z2: enzyme, enzyme_strand, enzyme_pos
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('-i',dest='i',help='input BAM file', required=True)
parser.add_argument('-o',dest='o',help='output BAM file', required=True)
parser.add_argument('-b',dest='b',help='tabix-indexed GBSBED file of restriction sites', required=True)
args = parser.parse_args()


def fetch_sites(chrom, bed):
    '''Fetch restriction sites and store in a dictionary.'''
    # sites[chrom] is keyed by (position, read.is_reverse).
    sites = {chrom:{}}
    try:
        for line in bed.fetch(chrom):
            line = line.strip()
            fields = line.split()
            try:
                pos = int(fields[1])
            except:
                pos = None
            enzyme, strand = fields[5:7]
            fwd_sites, rev_sites = fields[10:]
            fwd_insert, rev_insert = fields[3:5]
            try:
                fwd_insert = int(fwd_insert)
            except:
                fwd_insert = None
            try:
                rev_insert = int(rev_insert)
            except:
                rev_insert = None
            # forward read ligation sites (0-indexed).
            fwd_sites = [int(i) for i in fwd_sites.split(',')]
            # reverse read ligation sites.
            rev_sites = [int(i) for i in rev_sites.split(',')]
            for site in rev_sites:
                sites[chrom][(site, True)] = (enzyme, strand, pos, rev_insert)
            for site in fwd_sites:
                sites[chrom][(site, False)] = (enzyme, strand, pos, fwd_insert)
    except:
        pass
    return(sites)


# Bed file of restriction sites, indexed by tabix.
restrictbed = pysam.Tabixfile(args.b, 'r')
inbam = pysam.Samfile(args.i, 'rb')
outbam = pysam.Samfile(args.o, 'wb', template=inbam)
# Each read is queried against a hash of restriction sites.
restrict = {}
# Read counter.
n = 0

# Parse the bam file and check if each read maps to a restriction site.
for read in inbam.fetch():
    n += 1
    if n % 100000 == 0:
        print '%i reads processed' % n
        stdout.flush()
    ligation_site = None
    enzyme = None
    strand = None
    pos = None
    insert = None
    if not read.is_unmapped:
        chrom = inbam.getrname(read.tid)
        # Make hash of restriction sites.
        if not restrict or chrom not in restrict:
            restrict = fetch_sites(chrom, restrictbed)
        # Look up the ligation site in the hash.
        if read.is_reverse:
            # Account for dynamic trimming.
            ligation_site = (read.aend - 1) + (read.rlen - read.qend)
            try:
                enzyme, strand, pos, insert = restrict[chrom][(ligation_site, True)]
            except:
                pass
        else:
            # Account for dynamic trimming.
            ligation_site = read.pos - read.qstart
            try:
                enzyme, strand, pos, insert = restrict[chrom][(ligation_site, False)]
            except:
                pass
    # Update the read tags.
    if insert:
        read.tags += [('Z0', insert)]
    if ligation_site:
        read.tags += [('Z1', int(ligation_site))]
    if enzyme:
        enzyme_tag = "%s;%s;%i" % (enzyme, strand, pos)
        read.tags += [('Z2', enzyme_tag)]
    outbam.write(read)
print "%i reads processed" % n
print "done!"
inbam.close()
outbam.close()
