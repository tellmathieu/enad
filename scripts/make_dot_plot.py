from Bio import SeqIO
from Bio.SeqUtils import nt_search

'''def getArgs():
  length = sys.argv[1]

  fasta1 = sys.argv[1]
  fasta2 = sys.argv[2]
  fasta3 = sys.argv[3]
  fasta4
  return'''


fasta1 = "/group/ctbrowngrp/finnolab_shared/eNAD_PacBio/resources/reference_genome_chr_only.fa"
fasta2 = "/group/ctbrowngrp/finnolab_shared/eNAD_PacBio/horse_01_200/horse_01_200_scaffold_chr.fa"

genome1 = list(SeqIO.parse(fasta1, "fasta"))
genome2 = list(SeqIO.parse(fasta2, "fasta"))

record1 = genome1[0]
record2 = genome2[0]

sequence1 = str(record1.seq)
sequence2 = str(record2.seq)

# Set parameters
window_size = 10
threshold = 7

# Generate DOT plot
dot_plot = nt_search(sequence1, sequence2)

# Display or save the DOT plot
import matplotlib.pyplot as plt
plt.imshow(dot_plot, cmap='Greys', interpolation='none')
plt.title('DOT Plot')
plt.xlabel('Sequence 2 Position')
plt.ylabel('Sequence 1 Position')
plt.savefig('test.png')
