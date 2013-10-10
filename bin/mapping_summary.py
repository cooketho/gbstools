from optparse import OptionParser
import pysam

"""Summarize restriction site mapping for a bam file.
"""

# Parse command line arguments from user.
USAGE = """
mapping_summary.py -i <input bam file>
"""
parser = OptionParser(USAGE)
parser.add_option('-i',dest='i',help='input bam file')
(options,args)=parser.parse_args()

bam = pysam.Samfile(options.i, 'rb')
enzyme_counts = {}
insert_counts = {}
# Keyed by: no. mapped, no. mapped to restriction site.
read_counts = {0:{0:0, 1:0, 2:0},
               1:{0:0, 1:0, 2:0},
               2:{0:0, 1:0, 2:0}}

print "# file: %s" % options.i
print "#"
for read in bam.fetch():
    tags = dict(read.tags)
    enzyme1, enzyme2, insert = None, None, None
    try:
        # Get the enzyme tag.
        enzyme1 = tags['Z2'].split(';')[0]
        if enzyme1 in enzyme_counts:
            # Increment the read counter for that enzyme.
            enzyme_counts[enzyme1] += 1
        else:
            enzyme_counts[enzyme1] = 1
    except:
        pass
    try:
        # Get the mate enzyme tag.
        enzyme2 = tags['Z4'].split(';')[0]
    except:
        pass
    try:
        # Get the insert size.
        insert = tags['Z0']
        if insert in insert_counts:
            # Increment the read counter for that insert size.
            insert_counts[insert] += 1
        else:
            insert_counts[insert] = 1
    except:
        pass
    if read.is_paired:
        # Number of reads in pair that mapped.
        mapped = 2 - (read.is_unmapped + read.mate_is_unmapped)
    else:
        mapped = not read.is_unmapped
    # Number of read in pair that mapped to a restriction site.
    enzyme_mapped = bool(enzyme1) + bool(enzyme2)
    # Increment the read counter.
    read_counts[mapped][enzyme_mapped] += 1

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
