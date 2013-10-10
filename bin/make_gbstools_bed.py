import sys
from optparse import OptionParser
from collections import namedtuple

"""Convert a bed file of restriction sites into a bed-like format used
by annotate_reads.py and other GBStools scripts.
"""

# Parse command line arguments from user.
USAGE = """
make_gbstools_bed.py -i <input bed file>
                     -o <output GBSbed file> (default=stdout)
                     --sort_by_cutsite <bool> (default=False)
"""

parser = OptionParser(USAGE)
parser.add_option('-i', '--input', dest='i', default=sys.stdin, help='input bed file')
parser.add_option('-o', '--output', dest='o', default=sys.stdout, help='output GBSbed file')
parser.add_option('--sort_by_cutsite', dest='sort', action="store_true", help='sort by cut site?')
(opt, args) = parser.parse_args()

# Subtract these values from the restriction site start to get the expected ligation sites.
LIG_SITE = {True:{'+':{'BpuEI':[-20,-19],
                       'BsaXI':[-19,-18,13,14],
                       'CspCI':[-23,-22,13,14],
                       'ApeKI':[-3]},
                  '-':{'BpuEI':[16,17],
                       'BsaXI':[-21,-20,11,12],
                       'CspCI':[-23,-22,13,14]}},
            False:{'+':{'BpuEI':[-22,-21],
                        'BsaXI':[-21,-22,9,10],
                        'CspCI':[-25,-24,10,11],
                        'ApeKI':[-1]},
                   '-':{'BpuEI':[13,14],
                        'BsaXI':[-24,-23,7,8],
                        'CspCI':[-25,-24,10,11]}}}

RecSite = namedtuple('RecSite', 'enzyme start end strand')
chrom = None
header = ['chrom', 'start', 'end', 'fwd_insert', 'rev_insert', 'enzyme', 'strand',
          'id', 'freq', 'cut_allele', 'fwd_ligation_sites', 'rev_ligation_sites']
print '#' + '\t'.join(header)

class Chromosome():
    # Class for storing and manipulating restriction site data from one chromosome.
    def __init__(self, chrom):
        self.chrom = chrom
        self.sites = []
        self.sites_merged = None

    def merge_sites(self):
        # Sort by the max of the rev ligation sites or the min of the forward sites.
        if opt.sort:
            sites_sorted = sorted(self.sites, key=lambda i: i['frag_end'])
        else:
            sites_sorted = self.sites
        # Loop through sorted list and assign an expected GBS fragment length 
        # (only consider pairs where site i is fwd and site i+1 is rev).
        for i in range(len(sites_sorted) - 1):
            if 'fwd_lig' in sites_sorted[i] and 'rev_lig' in sites_sorted[i+1]:
                frag_len = sites_sorted[i+1]['frag_end'] - sites_sorted[i]['frag_end'] + 1
                sites_sorted[i]['fwd_len'] = frag_len
                sites_sorted[i+1]['rev_len'] = frag_len
        # Merge fwd and rev sites by enzyme.
        sites_merged = {}
        for site in sites_sorted:
            site_key = site['RecSite']
            if site_key not in sites_merged:
                # Hash keyed by RecSite.
                sites_merged[site_key] = {}
            if 'fwd_lig' in site:
                sites_merged[site_key]['fwd_lig'] = site['fwd_lig']
                sites_merged[site_key]['fwd_len'] = site['fwd_len']
            else:
                sites_merged[site_key]['rev_lig'] = site['rev_lig']
                sites_merged[site_key]['rev_len'] = site['rev_len']
        self.sites_merged = sites_merged
        return(None)

    def write(self):
        sites = sorted(self.sites_merged.keys(), key=lambda s: s.start)
        for site in sites:
            # Arrange output in a list.
            fwd_lig = ','.join([str(i) for i in self.sites_merged[site]['fwd_lig']])
            rev_lig = ','.join([str(i) for i in self.sites_merged[site]['rev_lig']])
            output = [self.chrom, site.start, site.end, 
                      self.sites_merged[site]['fwd_len'],
                      self.sites_merged[site]['rev_len'], 
                      site.enzyme, site.strand, '.', '.', 'ref', fwd_lig, rev_lig]
            print '\t'.join([str(i) for i in output])
        return(None)
      
  
if opt.i != sys.stdin:
    input = open(opt.i, 'r')
else:
    input = sys.stdin

chrom = None
# Build a list of expected ligation sites for fwd and rev reads for each enzyme site.
for line in input:
    line = line.strip()
    if line[0] == "#":
        continue
    # Unpack the fields in the digest file.
    chrom_name, start, end, enzyme, score, strand = line.split()
    start, end = int(start), int(end)
    if not chrom:
        chrom = Chromosome(chrom_name)
    if chrom_name != chrom.chrom:    # Print out the sites from the last chromosome.
        chrom.merge_sites()
        chrom.write()
        chrom = Chromosome(chrom_name)    # Create a new Chromosome object.
    # Subtract 1 from the start site to convert to 0-indexed start site.
    fwd_lig = [start - i for i in LIG_SITE[False][strand][enzyme]]
    rev_lig = [start - i for i in LIG_SITE[True][strand][enzyme]]
    # Append the site to the Chromosome object (fwd and rev).
    chrom.sites.append({'RecSite':RecSite(enzyme, start, end, strand),
                        'rev_lig':rev_lig, 'frag_end':min(rev_lig), 'rev_len':'.'})
    chrom.sites.append({'RecSite':RecSite(enzyme, start, end, strand),
                        'fwd_lig':fwd_lig, 'frag_end':max(fwd_lig), 'fwd_len':'.'})                  
# Print out the last chromosome in the file.
chrom.merge_sites()
chrom.write()
