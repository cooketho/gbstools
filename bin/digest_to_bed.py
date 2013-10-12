#!/usr/bin/env python
from sys import stdin

"""Convert emboss restriction output to BED format.

Usage: cat <myfile.digest> | python digest_to_bed.py > <myfile.bed>

If emboss restrict uses the name of a different enzyme than the
one you specified and you want to change it to the correct name, use sed
to make the correction, e.g. to change ''TseI'' to ''ApeKI'' use:

cat <myfile.digest> | sed 's/TseI/ApeKI/' | python digest_to_bed.py > <myfile.bed>
"""

header = ('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand')
print '#' + '\t'.join(header)

for line in stdin:
    line = line.strip()
    fields = line.split()
    # Ignore blank lines, comments, and the header.
    if not line or line[:5]=="Start":
        continue
    elif line[0]=="#":
        if "Sequence:" in fields and "from:" in fields and "to:" in fields:
            chrom = fields[2]
    else:
        start, end, strand, enzyme = fields[:4]
        start = int(start)
        output = '\t'.join((chrom, str(start - 1), end, enzyme, '.', strand))
        print output
