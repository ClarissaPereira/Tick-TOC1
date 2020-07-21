# searches for motifs that appear in each separate DNA string (taking into account up to d mutations) 

import itertools

def all_kmers(k):
    return (''.join(p) for p in itertools.product('ATCG', repeat=k))

def hamming_distance(pattern, seq):
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(pattern, seq))
# zipping the strings to be compared into tuples cuts down execution time

def sub_strings(string, k):
    for i in range(len(string) - k + 1):
        yield string[i:i+k]
# making this a function will eliminate need for nested for loops

def basic_motif_search(k, d, DNA):
    motifs = set()
    for pos_kmer in all_kmers(k):
        if all(any(hamming_distance(pos_kmer, sub_string) <= d for sub_string in sub_strings(string, k)) for string in DNA):
# using the functions 'any' and 'all' will short circuit the execution and so cut down execution time - i.e. won't keep calculating hamming distance if it exceeds d 
# flattened the previous nested loop
            motifs.add(pos_kmer)
    return motifs

DNA = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']

print(basic_motif_search(3, 1,DNA))
