# mapping_summary.py
# Tom Cooke
# 2013-08-22
# Summarize mapping of reads to restriction sites. See README for details about format.

from optparse import OptionParser
import pysam

# Parse command line arguments from user.
USAGE = """
rs_coverage_histogram.py -i <input bam file>
"""
parser = OptionParser(USAGE)
parser.add_option('-i',dest='i',help='input bam file')
(options,args)=parser.parse_args()

bam = pysam.Samfile(options.i, 'rb')
enzyme_counts = {}
insert_counts = {}
read_counts = {0:{0:0, 1:0, 2:0},
               1:{0:0, 1:0, 2:0},
               2:{0:0, 1:0, 2:0}}    # Keyed by: no. mapped, no. mapped to restriction site.

print "# file: %s" % options.i
print "#"
for read in bam.fetch():
    tags = dict(read.tags)
    enzyme1, enzyme2, insert = None, None, None
    try:
        enzyme1 = tags['Z2'].split(';')[0]    # Get the enzyme tag.
        if enzyme1 in enzyme_counts:
            enzyme_counts[enzyme1] += 1    # Increment the read counter for that enzyme.
        else:
            enzyme_counts[enzyme1] = 1
    except:
        pass
    try:
        enzyme2 = tags['Z4'].split(';')[0]    # Get the mate enzyme tag.
    except:
        pass
    try:
        insert = tags['Z0']    # Get the insert size.
        if insert in insert_counts:
            insert_counts[insert] += 1    # Increment the read counter for that insert size.
        else:
            insert_counts[insert] = 1
    except:
        pass
    mapped = 2 - (read.is_unmapped + read.mate_is_unmapped)    # Number of reads in pair that mapped.
    enzyme_mapped = bool(enzyme1) + bool(enzyme2)    # Number of read in pair that mapped to a restriction site.
    read_counts[mapped][enzyme_mapped] += 1    # Increment the read counter.

print "# mapped_in_pair\tmapped_in_pair_to_restriction_site\treads"
for i in (0, 1, 2):
    for j in (0, 1, 2):
        print "# %i\t%i\t%i" % (i, j, read_counts[i][j])
print "#"

print "# enzyme\treads"
for enzyme in enzyme_counts:
    print "# %s\t%i" % (enzyme, enzyme_counts[enzyme])
print "#"

print "insert_size\treads"
inserts = sorted(insert_counts.keys())
for insert in inserts:
    print "%i\t%i" % (insert, insert_counts[insert])
