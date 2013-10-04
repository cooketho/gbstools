# annotate_bam.py
# Tom Cooke
# 2013-08-20
# For every read in a bam file add the following tags:
# Z0: ligation site
# Z1: enzyme, enzyme_strand, enzyme_pos
# Z2: mate ligation site
# Z3: mate_enzyme, mate_enzyme_strand, mate_enzyme_pos

import pysam
from sys import stdout
from optparse import OptionParser
from collections import namedtuple

# Function to fetch restriction sites from one chromosome and store in a hash
def fetch_sites(chrom, bed):
    sites = {chrom:{True:{}, False:{}}}    # Key by chrom, then "read.is_reverse" boolean, then by position.
    try:
        for line in bed.fetch(chrom):
            line = line.strip()
            fields = line.split()
            pos = fields[1]
            enzyme, strand = fields[5:7]
            fwd_sites, rev_sites = fields[10:]
            fwd_sites = [int(i) for i in fwd_sites.split(',')]    # forward read ligation sites (0-indexed).
            rev_sites = [int(i) for i in rev_sites.split(',')]    # reverse read ligation sites.
            for site in rev_sites:
                if site not in sites[chrom][True]:
                    sites[chrom][True][site] = (enzyme, strand, pos)
            for site in fwd_sites:
                if site not in sites[chrom][False]:
                    sites[chrom][False][site] = (enzyme, strand, pos)
    except:
        pass
    return(sites)

# Parse command line arguments from user.
USAGE = """
annotate_bam.py -i <input bam file>
                -o <output bam file>
                -r <restriction site bed file>
format for bed file (tab-separated with comma-separated sub-fields):
chrom, start, end, fwd_insert, rev_insert, enzyme, strand, pos, id, freq, cut_allele, fwd_ligation_sites, rev_ligation_sites
"""
parser = OptionParser(USAGE)
parser.add_option('-i',dest='i',help='input bam file')
parser.add_option('-o',dest='o',help='output bam file')
parser.add_option('-r',dest='r',help='restriction site bed file')
(options,args)=parser.parse_args()

restrictbed = pysam.Tabixfile(options.r, 'r')    # Bed file of restriction sites, indexed by tabix.
inbam = pysam.Samfile(options.i, 'rb')
outbam = pysam.Samfile(options.o, 'wb', template=inbam)
stack = []    # Reads are put in a stack until the mate pair is found, or the stack becomes too big.
stack_range = 2000    # Pop reads out of stack if the stack range is > 2000 bp.
reads = {}    # Hash for fast lookup of mate pairs in stack.
restrict = {}    # Each read is queried against a hash of restriction sites.
n = 0    # Read counter.

Read = namedtuple('Read', 'read tags ligation_site enzyme_tag')

# Parse the bam file and check if each read maps to a restriction site.
for read in inbam.fetch('1', 1, 1000000):
    n += 1
    if n % 100000 == 0:
        print '%i reads processed' % n
        stdout.flush()
    qname = read.qname.split('_')[0]    # In the newer Illumina format the mate pair identifier is the part before '_'
    is_read1 = read.is_read1    # Is this the first or second read in the mate pair?
    ligation_site, enzyme, strand, pos = [None] * 4    # Defaults.
    tags = []
    if not read.is_unmapped:
        chrom = inbam.getrname(read.tid)
        # If necessary, make a new hash of restriction sites on this chrom, keyed by: chrom, read.is_reverse, ligation site.
        if not restrict or chrom not in restrict:
            restrict = fetch_sites(chrom, restrictbed)
        # Check if the ligation site is in the hash.
        if read.is_reverse:
            ligation_site = (read.aend - 1) + (read.rlen - read.qend)    # Account for dynamic trimming.
            if ligation_site in restrict[chrom][True]:
                enzyme, strand, pos = restrict[chrom][True][ligation_site]
        else:
            ligation_site = read.pos - read.qstart    # Account for dynamic trimming.
            if ligation_site in restrict[chrom][False]:
                enzyme, strand, pos = restrict[chrom][False][ligation_site]
    # Update the read tags.
    if ligation_site:
        tags.append(('Z1', int(ligation_site)))
    if enzyme:
        enzyme_tag = ';'.join([str(i) for i in (enzyme, strand, pos)])
        tags.append(('Z2', enzyme_tag))
    else:
        enzyme_tag = None

    # For reverse reads, check if the mate pair is in the stack.
    if read.is_reverse and (read.qname, read.is_read2) in reads:
        mate = reads[(read.qname, read.is_read2)]
        # Add the mate pair ligation site to the tags if it exists.
        if mate.ligation_site:
            tags.append(('Z3', int(mate.ligation_site)))    # Convert from long to int.
            if ligation_site:
                # Calculate the insert size.
                insert = abs(mate.ligation_site - ligation_site)
                tags.append(('Z0', int(insert)))    # Convert from long to int.
        # Add the mate pair enzyme tag if it exists.
        if mate.enzyme_tag:
            tags.append(('Z4', mate.enzyme_tag))

    # Add the read to the stack.
    stack.append((qname, read.is_read1))
    reads[(qname, read.is_read1)] = Read(read, tags, ligation_site, enzyme_tag)

    # If the bottom read in the stack is unmapped or > stack_range bp away from the top, remove it.
    br = reads[stack[0]]
    if (br.read.is_unmapped or 
        (not read.is_unmapped and (br.read.tid != read.tid or br.read.pos < read.pos - stack_range))):
        # For forward reads, check if the mate pair is in the stack.
        if not br.read.is_reverse and (br.read.qname, br.read.is_read2) in reads:
            mate = reads[(br.read.qname, br.read.is_read2)]
            # Add the mate pair ligation site to the tags if it exists.
            if mate.ligation_site:
                br.tags.append(('Z3', int(mate.ligation_site)))
                # Calculate the insert size.
                if br.ligation_site:
                    insert = abs(mate.ligation_site - br.ligation_site)
                    br.tags.append(('Z0', int(insert)))
            # Add the mate pair enzyme tag if it exists.
            if mate.enzyme_tag:
                br.tags.append(('Z4', mate.enzyme_tag))
        br.read.tags += sorted(br.tags)
        # Write the bottom read to file and pop it from the stack.
        outbam.write(br.read)
        del reads[stack[0]]
        del stack[0]

# After parsing the whole file, empty the stack.
for readname in stack:
    br = reads[readname]    # Get the read object itself.
    # Check if the mate pair is in the stack.
    if not br.read.is_reverse and (br.read.qname, br.read.is_read2) in reads:
        mate = reads[(br.read.qname, br.read.is_read2)]
        # Add the mate pair ligation site to the tags if it exists.
        if mate.ligation_site:
            br.tags.append(('Z3', int(mate.ligation_site)))
            # Calculate the insert size.
            if br.ligation_site:
                insert = abs(mate.ligation_site - br.ligation_site)
                br.tags.append(('Z0', int(insert)))
        # Add the mate pair enzyme tag if it exists.
        if mate.enzyme_tag:
            br.tags.append(('Z4', mate.enzyme_tag))
    br.read.tags += sorted(br.tags)
    outbam.write(br.read)
print "%i reads processed" % n
print "done!"
inbam.close()
outbam.close()
