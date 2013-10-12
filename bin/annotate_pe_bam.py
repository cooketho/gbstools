#!/usr/bin/env python
import pysam
from sys import stdout
import argparse
from collections import namedtuple

# Parse command line arguments from user.
USAGE = """
annotate_pe_bam.py -i <input bam file>
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
Z3: mate ligation site (mapping position of the opposite end of the SAM ``Template``)
Z4: mate_enzyme, mate_enzyme_strand, mate_enzyme_pos
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('-i',dest='i',help='input BAM file', required=True)
parser.add_argument('-o',dest='o',help='output BAM file', required=True)
parser.add_argument('-b',dest='b',help='tabix-indexed GBSBED file of restriction sites', required=True)
args = parser.parse_args()


def update_tags(stack_read, reads):
    '''Update a StackRead with mate pair info from the stack of reads.'''
    try:
        # Mate has same qname but opposite ''is_read1'' bool.
        mate = reads[(stack_read.read.qname, not stack_read.read.is_read1)]
        try:
            stack_read.tags.append(('Z3', int(mate.ligation_site)))
        except:
            pass
        try:
            insert = abs(mate.ligation_site - stack_read.ligation_site)
            stack_read.tags.append(('Z0', int(insert)))
        except:
            pass
        try:
            stack_read.tags.append(('Z4', mate.enzyme_tag))
        except:
            pass
    except:
        pass
    return(None)


def fetch_sites(chrom, bed):
    '''Fetch restriction sites and store in a dictionary.'''
    # sites[chrom] is keyed by (position, read.is_reverse).
    sites = {chrom:{}}
    try:
        for line in bed.fetch(chrom):
            line = line.strip()
            fields = line.split()
            pos = int(fields[1])
            enzyme, strand = fields[5:7]
            fwd_sites, rev_sites = fields[10:]
            # forward read ligation sites (0-indexed).
            fwd_sites = [int(i) for i in fwd_sites.split(',')]
            # reverse read ligation sites.
            rev_sites = [int(i) for i in rev_sites.split(',')]
            for site in rev_sites:
                sites[chrom][(site, True)] = (enzyme, strand, pos)
            for site in fwd_sites:
                sites[chrom][(site, False)] = (enzyme, strand, pos)
    except:
        pass
    return(sites)


# Bed file of restriction sites, indexed by tabix.
restrictbed = pysam.Tabixfile(args.b, 'r')
inbam = pysam.Samfile(args.i, 'rb')
outbam = pysam.Samfile(args.o, 'wb', template=inbam)
# Stack of reads is used to speed up mate pair searches.
stack = []
# Max stack size.
stack_range = 2000
# Hash for fast lookup of mate pairs in stack.
reads = {}
# Each read is queried against a hash of restriction sites.
restrict = {}
# Read counter.
n = 0

StackRead = namedtuple('StackRead', 'read tags ligation_site enzyme_tag')

# Parse the bam file and check if each read maps to a restriction site.
for read in inbam.fetch():
    n += 1
    if n % 100000 == 0:
        print '%i reads processed' % n
        stdout.flush()
    # For Cassava >= v1.8 the " " separator is converted to "_" by bwa.
    qname = read.qname.split('_')[0]
    # Is this read 1 or 2 in the mate pair?
    is_read1 = read.is_read1    
    ligation_site = None
    enzyme = None
    strand = None
    pos = None
    tags = []
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
                enzyme, strand, pos = restrict[chrom][(ligation_site, True)]
            except:
                pass
        else:
            # Account for dynamic trimming.
            ligation_site = read.pos - read.qstart
            try:
                enzyme, strand, pos = restrict[chrom][(ligation_site, False)]
            except:
                pass
    # Update the read tags.
    if ligation_site:
        tags.append(('Z1', int(ligation_site)))
    if enzyme:
        enzyme_tag = "%s;%s;%i" % (enzyme, strand, pos)
        tags.append(('Z2', enzyme_tag))
    else:
        enzyme_tag = None

    # Create a ''StackRead'' to store the read data.
    stack_read = StackRead(read, tags, ligation_site, enzyme_tag)
    # Look up mate pair tags.
    if stack_read.read.is_reverse:
        update_tags(stack_read, reads)

    # Add the read to the stack.
    stack.append((qname, read.is_read1))
    reads[(qname, read.is_read1)] = stack_read

    # Get the top and bottom reads from the stack.
    top = reads[stack[-1]]
    bottom = reads[stack[0]]
    # Should the bottom read be popped from the stack?
    if bottom.read.is_unmapped:
        pop = True
    elif bottom.read.tid != top.read.tid:
        pop = True
    elif bottom.read.pos < top.read.pos - stack_range:
        pop = True
    else:
        pop = False
    if pop:
        # If the bottom read is forward-mapped, look for the mate in the stack.
        if not bottom.read.is_reverse:
            update_tags(bottom, reads)
        # Append the StackRead tags to the pysam read object.
        bottom.read.tags += sorted(bottom.tags)
        # Write the read and then remove it from the stack.
        outbam.write(bottom.read)
        del reads[stack[0]]
        del stack[0]
# When done parsing the input, empty the stack.
for read in stack:
    stack_read = reads[read]
    update_tags(stack_read, reads)
    stack_read.read.tags += sorted(stack_read.tags)
    outbam.write(stack_read.read)
print "%i reads processed" % n
print "done!"
inbam.close()
outbam.close()
