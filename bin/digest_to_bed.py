# digest_to_bed.py
# Tom Cooke
# 2013-09-04
# Convert emboss restrict output to bed-like file format for GBStools.

from sys import stdin
from collections import namedtuple

# Subtract these values from the restriction site start to get the expected ligation sites.
expected_sites = {True:{'+':{'BpuEI':[-20,-19],
                            'BsaXI':[-19,-18,13,14],
                            'CspCI':[-23,-22,13,14]},
                       '-':{'BpuEI':[16,17],
                            'BsaXI':[-21,-20,11,12],
                            'CspCI':[-23,-22,13,14]}},
                 False:{'+':{'BpuEI':[-22,-21],
                             'BsaXI':[-21,-22,9,10],
                             'CspCI':[-25,-24,10,11]},
                        '-':{'BpuEI':[13,14],
                             'BsaXI':[-24,-23,7,8],
                             'CspCI':[-25,-24,10,11]}}}

RecognitionSite = namedtuple('RecognitionSite', 'enzyme start end strand')
chrom = None
header = ['chrom', 'start', 'end', 'fwd_insert', 'rev_insert', 'enzyme', 'strand',
          'id', 'freq', 'cut_allele', 'fwd_ligation_sites', 'rev_ligation_sites']
print '#' + '\t'.join(header)

class Chromosome():
    # Class for storing and manipulating restriction site data from one chromosome.
    def __init__(self, chrom):
        self.chrom = chrom
        self.sites = []

    def merge_sites(self):
        # Sort by the maximum of the rev ligation sites or the min of the forward sites.
        sites_sorted = sorted(self.sites, key=lambda i: i['sortkey'])
        # Loop through sorted list and assign an expected GBS fragment length 
        # (only consider pairs where site i is fwd and site i+1 is rev).
        # If this is not the case the length defaults to '.'.
        for i in range(len(sites_sorted)-1):
            if 'fwd' in sites_sorted[i] and 'rev' in sites_sorted[i+1]:
                frag_len = sites_sorted[i+1]['sortkey'] - sites_sorted[i]['sortkey']
            if 'fwd' in sites_sorted[i]:
                sites_sorted[i]['fwd_len'] = frag_len
            if 'rev' in sites_sorted[i+1]:
                sites_sorted[i+1]['rev_len'] = frag_len
        # Merge fwd and rev sites by enzyme.
        sites_merged = {}
        for site in sites_sorted:
            site_key = site['recognition_site']
            if site_key not in sites_merged:
                sites_merged[site_key] = {}    # Hash keyed by recognition site.
            if 'fwd' in site:
                sites_merged[site_key]['fwd'] = site['fwd']    # Forward read recognition sites_merged.
                sites_merged[site_key]['fwd_len'] = site['fwd_len']    # Forward read insert length.
            else:
                sites_merged[site_key]['rev'] = site['rev']    # Reverse read recognition sites_merged.
                sites_merged[site_key]['rev_len'] = site['rev_len']    # Reverse read insert length.
        # Sort keys by recognition site.
        self.recog_sites = sorted(sites_merged.keys(), key=lambda i: i.start)
        self.sites_merged = sites_merged
        return(None)

    def write(self):
        # Print one line per enzyme in a format that the count_reads.py script can understand.
        for recog_site in self.recog_sites:
            # Arrange output in a list.
            fwd_ligation_sites = ','.join([str(i) for i in self.sites_merged[recog_site]['fwd']])
            rev_ligation_sites = ','.join([str(i) for i in self.sites_merged[recog_site]['rev']])
            output = [self.chrom, recog_site.start, recog_site.end, self.sites_merged[recog_site]['fwd_len'],
                      self.sites_merged[recog_site]['rev_len'], recog_site.enzyme, recog_site.strand, 
                      '.', '.', 'ref', fwd_ligation_sites, rev_ligation_sites]
            print '\t'.join([str(i) for i in output])
        return(None)
        

# Build a list of expected ligation sites for fwd and rev reads for each enzyme site.
for line in stdin:
    line = line.strip()
    # Ignore blank lines, comments, and the header.
    if line and line[0]=="#":
        fields = line.split()
        if "Sequence:" in fields and "from:" in fields and "to:" in fields:
            if chrom:    # Print out the sites from the last chromosome.
                chrom.merge_sites()
                chrom.write()
            chrom = Chromosome(fields[2])    # Create a new Chromosome object.
        continue
    elif len(line)==0 or line[:5]=="Start":
        continue
    # Unpack the fields in the digest file.
    start, end, strand, enzyme, sequence, five_prime, three_prime, five_primeR, three_primeR = line.split()
    start, end = int(start), int(end)
    # Subtract 1 from the start site to convert to 0-indexed start site.
    fwd_ligation_sites = [(start - 1) - i for i in expected_sites[False][strand][enzyme]]
    rev_ligation_sites = [(start - 1) - i for i in expected_sites[True][strand][enzyme]]
    # Append the site to the Chromosome object (fwd and rev).
    chrom.sites.append({'recognition_site':RecognitionSite(enzyme, start - 1, end, strand),
                  'fwd':fwd_ligation_sites, 'sortkey':max(fwd_ligation_sites), 'fwd_len':'.'})                  
    chrom.sites.append({'recognition_site':RecognitionSite(enzyme, start - 1, end, strand),
                  'rev':rev_ligation_sites, 'sortkey':min(rev_ligation_sites), 'rev_len':'.'})
# Print out the last chromosome in the file.
chrom.merge_sites()
chrom.write()
