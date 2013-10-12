#!/usr/bin/env python
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
import argparse

USAGE = """
simulate_fasta.py -N <sequence length> > <myfile.fa>
"""

DESCRIPTION = """
Draw a random sequence of length N and output to stdout in fasta format.
(requires BioPython).
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('-N', dest='N', type=int, help='number of bases', required=True)
parser.add_argument('--seed', dest='seed', type=int, default=0, help='seed for random number generator')
parser.add_argument('--name', dest='name', type=str, default='S', help='sequence name')
parser.add_argument('--description', dest='description', type=str, default='simulated_chromosome', help='sequence description')
args = parser.parse_args()

N = args.N
rand = random.Random(args.seed)
ref = [rand.choice(('A', 'G', 'T', 'C')) for i in range(N)]
ref = Seq(''.join(ref), IUPAC.unambiguous_dna)
refseq = SeqRecord(ref, id=args.name, description=args.description)
SeqIO.write([refseq], sys.stdout, 'fasta')
