from optparse import OptionParser

"""Parse RS mapping summaries from mapping_summary.py and create
a table of depth-of-coverage normalization factors to be used
with GBStools.
"""

# Parse command line arguments from user.
USAGE = """
normfactors.py --summaries <summarylist> > nf.txt
(summary list is in sample/summary_file format)
"""

parser = OptionParser(USAGE)
parser.add_option('--summaries', dest='summaries', help='list of mapping summary files')
parser.add_option('--max_insert', dest='max_insert', type="int", default=1000, help='maximum insert size')
parser.add_option('--window', dest='window', type="int", default=0, help='sliding window size for smoothing')
(opt, args) = parser.parse_args()

max_insert = opt.max_insert
win = opt.window

# Parse the rs mapping summaries and store insert counts in a dict.
counts = {}
summaries = open(opt.summaries, 'r')
for line in summaries:
    line = line.strip()
    sample, summary = line.split()
    # Store counts in a list.
    counts[sample] = [0] * max_insert
    summary = open(summary, 'r')
    for line in summary:
        line = line.strip()
        # Ignore headers.
        if not line or line[0] == "#":
            continue
        line = line.split()
        if line[0] == "insert_size":
            continue
        insert = int(line[0])
        count = int(line[1])
        try:
            counts[sample][insert] = count
        except:
            pass

# Normalize counts for each insert size.
samples = sorted(counts.keys())
normfactors = []
for insert in range(max_insert):
    count = [counts[sample][insert] for sample in samples]
    mean_count = float(sum(count)) / len(count)
    try:
        nf = [c / mean_count for c in count]
    except:
        nf = [0] * len(count)
    normfactors.append(dict(zip(samples, nf)))

# Print header line.
print "#insert\t" + "\t".join(samples)
# Take a sliding-window mean of the normalization factors.
for insert in range(1, max_insert):
    window = range(max(insert - win, 0), min(insert + win + 1, max_insert))
    nf = {}
    for sample in samples:
        nf[sample] = sum([normfactors[i][sample] for i in window]) / len(window)
    output = "%i\t" % insert
    output = output + "\t".join(['{0:.3f}'.format(nf[sample]) for sample in samples])
    print output
