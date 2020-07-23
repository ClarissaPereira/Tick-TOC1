# This function evaluates the probability of every possible k-mer in a sequence and identifies the most probable k-mer
# Takes inputs of a position weight matrix (referred to as probability_matrix), a sequence of DNA, and k 

def generate_probable_kmers (probability_matrix, k):
  base_value = {0:'A', 1:'C', 2:'G', 3:'T'}
# links column number to nucleotide
  max_prob_dict = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[], 9:[], 10:[], 11:[], 12:[]}
# number of keys in dictionary must == k (this example is for a 13-mer so 13 keys in dictionary)
  most_prob_list = []
  for k in range(k):
      column = [column[k] for column in probability_matrix]
# creates a list for a column from the matrix
      maxc=(max(column))
      for i in range (len(column)):
          if column[i] == maxc:
              max_prob_dict[k].append(base_value[i]) 
  for k in range(k):
    most_prob_list.append(max_prob_dict[k])
  import itertools
  most_prob_kmers = [''.join(s) for s in itertools.product(*most_prob_list)]
# combines a list of lists of possible nucleotide combinations into a list of strings of probable k-mers
  return most_prob_kmers

def hamming_distance(most_prob_kmer, kmer):
  return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(most_prob_kmer, kmer))

def text_kmer_generator (sequence, k):
  text_kmers = []
  for n in range(len(sequence) - k + 1):
    text_kmers.append(sequence[n:n+k])
  return text_kmers

def most_probable_in_seq (sequence, k, probability_matrix): 
  most_prob_kmers = generate_probable_kmers(probability_matrix, k)
  text_kmers =  text_kmers(sequence, k)
  min_ham_dist = 100000
# pseudocode refers to infinity - using a very large number in place of infinity to avoid importing numpy                  
  for kmer in text_kmers:    
      for x in most_prob_kmers:
          if hamming_distance(most_prob_kmer, kmer) < min_ham_dist:
              min_ham_dist = hamming_distance(most_prob_kmer, kmer)
              best_match = kmer
print(best_match)

probability_matrix = [[0.2,0.2,0.3,0.2,0.3], [0.4,0.3,0.1,0.5,0.1], [0.3,0.3,0.5,0.2,0.4], [0.1,0.2,0.1,0.1,0.2]]
k = 5
sequence = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
most_probable_in_seq (sequence, k, probability_matrix)
# Output -> 'CCGAG'

