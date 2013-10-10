from optparse import OptionParser
import vcf
from Bio import SeqIO
import math
import random

"""Parse a VCF file of GBS SNPs, a bed file of restriction fragments
and a reference genome fasta and return a fastq file of reads for
each sample in the VCF.
"""

USAGE = """
make_gbstools_bed.py --vcf <snps vcf file>
                     --bed <GBStools restriction site bed file>
                     --fasta <reference genome fasta>
"""

parser = OptionParser(USAGE)
parser.add_option('--vcf', dest='vcf', help='snps vcf file')
parser.add_option('--bed', dest='bed', help='GBStools restriction site bed file')
parser.add_option('--fasta', dest='fasta', help='reference genomoe fasta')
(opt, args) = parser.parse_args()

rand = random.Random(0)
readlen = 101
maxlen = 500
epsilon = 0.001
def write_reads(seq, counts, myfile):
    '''Write sequence to fastq file.'''
    for i in range(counts):
        lane = rand.randint(1,8)
        tile = rand.randint(1,100)
        x = rand.randint(1,9999)
        y = rand.randint(1,9999)
        seqid = "@SIMULATION:%i:%i:%i:%i#0" % (lane, tile, x, y)
        qual = [chr(int(-10 * math.log(epsilon, 10) + 33))] * len(seq)
        qual = ''.join(qual)
        myfile.write(seqid + '\n')
        myfile.write(seq + '\n')
        myfile.write('+\n')
        myfile.write(qual + '\n')
        
# Get the reference sequence.
seq_record = SeqIO.parse(open(opt.fasta, 'r'), 'fasta').next()

# Get a list of reads from the in silico digest of the reference.
reads = []
bed = open(opt.bed, 'r')
for line in bed:
    line = line.strip()
    if line[0] == "#":
        continue
    (chrom, start, end, fwd_insert, rev_insert, enzyme, strand, 
     site_id, freq, cut_allele, fwd_lig, rev_lig) = line.split()
    if fwd_insert != '.' and rev_insert != '.':
        fwd_insert = int(fwd_insert)
        rev_insert = int(rev_insert)
        fwd_lig = max([int(i) for i in fwd_lig.split(',')])
        rev_lig = min([int(i) for i in rev_lig.split(',')])
        if fwd_insert >= 2 * readlen and fwd_insert <= maxlen:
            read = seq_record.seq[fwd_lig:fwd_lig + readlen]
            reads.append(str(read))
        if rev_insert >= 2 * readlen and rev_insert <= maxlen:
            read = seq_record.seq[rev_lig - readlen + 1:rev_lig + 1]
            read = read.reverse_complement()
            reads.append(str(read))
    else:
        pass
readiter = (read for read in reads)

# Make a hash of fastq files, keyed by sample name.
fastq = {}
reader = vcf.Reader(open(opt.vcf, 'r'))
for sample in reader.samples:
    fastq[sample] = open(sample + '.fastq', 'w')

# Parse the vcf file and write reads for each sample.
for snp in reader:
    ref = readiter.next()
    alt = [base for base in ref]
    bases = ['A', 'G', 'T', 'C']
    bases.remove(ref[50])
    alt[50] = rand.choice(bases)
    alt = ''.join(alt)
    for sample in snp.samples:
        if sample['AD']:
            write_reads(ref, int(sample['AD'][0]), fastq[sample.sample])
            write_reads(alt, int(sample['AD'][1]), fastq[sample.sample])
        else:
            pass

for sample in fastq:
    fastq[sample].close()
    
