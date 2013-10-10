import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from sys import argv, stdout

"""Draw a random sequence of length N and output to stdout in fasta format.

Usage: python simulate_fa.py <N> > <myfile.fa>

Requires BioPython
"""

script, N = argv

N = int(N)
rand = random.Random(0)
ref = [rand.choice(('A', 'G', 'T', 'C')) for i in range(N)]
ref = Seq(''.join(ref), IUPAC.unambiguous_dna)
refseq = SeqRecord(ref, id='S', description='simulated_chromosome')
SeqIO.write([refseq], stdout, 'fasta')
