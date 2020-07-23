probability_matrix = [[0.2,0.2,0.3,0.2,0.3], [0.4,0.3,0.1,0.5,0.1], [0.3,0.3,0.5,0.2,0.4], [0.1,0.2,0.1,0.1,0.2]
k = 8
base_value = {0:'A', 1:'C', 2:'G', 3:'T'}

max_prob_dict = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[], 9:[], 10:[], 11:[], 12:[], 13:[]}
most_prob_list = []
                      
for k in range(8):
    column = [column[k] for column in probability_matrix]
    maxc=(max(column))
    for i in range (len(column)):
        if column[i] == maxc:
            max_prob_dict[k].append(base_value[i])
print(max_prob_dict)

for k in range(8):
    most_prob_list.append(max_prob_dict[k])
import itertools
most_prob_kmers = [''.join(s) for s in itertools.product(*most_prob_list)]                     

def hamming_distance(most_prob_kmer, kmer):
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(most_prob_kmer, kmer))

sequence = '' 
min_ham_dist = 100000
text_kmers = []
for n in range(len(sequence) - 8 + 1):
    text_kmers.append(sequence[n:n+13]) 
                    
for kmer in text_kmers:    
    for x in most_prob_kmers:
        if hamming_distance(most_prob_kmer, kmer) < min_ham_dist:
            min_ham_dist = hamming_distance(most_prob_kmer, kmer)
            best_match = kmer
print(best_match)
