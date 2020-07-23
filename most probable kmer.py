probability_matrix = [[0.2,0.2,0.3,0.2,0.3], [0.4,0.3,0.1,0.5,0.1], [0.3,0.3,0.5,0.2,0.4], [0.1,0.2,0.1,0.1,0.2]
k = 8
base_value = {0:'A', 1:'C', 2:'G', 3:'T'}

most_prob_list = []
for i in range (k):
    motifp = 0
    for j in range (4):
        if (probability_matrix[j][i]) > motifp:
            motifp = (probability_matrix[j][i])
            base = base_value[j]
    most_prob_list.append(base)         
most_prob_kmer = (''.join(most_prob_list))
print (most_prob_kmer)

def hamming_distance(most_prob_kmer, kmer):
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(most_prob_kmer, kmer))

sequence = '' 
min_ham_dist = 100000
for n in range(len(sequence) - k + 1):
    kmer = sequence[n:n+k]
    if hamming_distance(most_prob_kmer, kmer) < min_ham_dist:
        min_ham_dist = hamming_distance(most_prob_kmer, kmer)
        print(min_ham_dist)
        best_match = kmer
print(best_match)
