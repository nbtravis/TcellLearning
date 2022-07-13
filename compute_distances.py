import pickle
import os

WORKING_DIR = os.path.dirname(os.path.realpath(__file__))
SELF_9MERS_FILE = os.path.join(WORKING_DIR, 'data/self-9mers.txt')

unique_self_peptides = set()
unique_aminos = set()
with open(SELF_9MERS_FILE, 'rb') as f:
    for line in f.readlines():
        peptide = line.decode('utf8').strip()
        unique_self_peptides.add(peptide)
        for amino in peptide:
            unique_aminos.add(amino)        

def peptides_with_one_edit(peptide):
    edited_peptides = set()
    for i in range(len(peptide)):
        amino = peptide[i]
        for other_amino in unique_aminos:
            if other_amino == amino:
                continue
            peptide_copy = list(peptide)
            peptide_copy[i] = other_amino
            edited_peptides.add(''.join(peptide_copy))
    return edited_peptides

def peptides_with_two_edits(peptide):
    edited_peptides = set()
    for i in range(len(peptide)):
        amino_i = peptide[i]
        for j in range(i + 1, len(peptide)):
            amino_j = peptide[j]

            for other_amino_i in unique_aminos:
                if other_amino_i == amino_i:
                    continue
                for other_amino_j in unique_aminos:
                    if other_amino_j == amino_j:
                        continue
                    peptide_copy = list(peptide)
                    peptide_copy[i] = other_amino_i
                    peptide_copy[j] = other_amino_j
                    edited_peptides.add(''.join(peptide_copy))
    return edited_peptides

# Create some sets of "edited" peptides for being able to quickly
# say what the nearest neighbor hamming distance is.
one_edit = set()
two_edits = set()

MODULUS=0
N = len(unique_self_peptides)
for i, self_peptide in enumerate(unique_self_peptides):
    if i % 10000 == 0:
        print('processing peptide {} / {}'.format(i, N))
    if i % 50 != MODULUS:
        continue
    one_edit |= peptides_with_one_edit(self_peptide)
    two_edits |= peptides_with_two_edits(self_peptide)


