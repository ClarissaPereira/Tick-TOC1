# Finds a k-mer motif which has the minimum hamming distance between itself and all given DNA sequences 

import itertools
def hamming_distance(pattern, seq):
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(pattern, seq))

def total_minimum_distances(DNA, pattern, k):
    total_distance = 0
    for text in DNA:
        distances = []
        sub_string = []
        for i in range(len(text)):
            sub_string.append(text[i:i+k])
        for j in sub_string:
            distances.append(hamming_distance(pattern, j))
        total_distance = total_distance + min(distances)
    return total_distance  

def all_kmers(k):
    return (''.join(p) for p in itertools.product('ATCG', repeat=k))

def median_motif(DNA, k):
    median_string = 'XXX'
    median = 10000
    for pos_kmer in all_kmers(k):
        if total_minimum_distances(DNA, pos_kmer, k) < median:
            median = total_minimum_distances(DNA, pos_kmer, k)
            median_string = pos_kmer
    return (median_string)
        
          
    
DNA = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG']

print(median_motif(DNA, 3))

# Output -> 'GAC'
