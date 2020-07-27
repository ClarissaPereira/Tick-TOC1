# generates n number of random k-mers from each DNA sequence in a list 

def random_kmers(DNA, t, k, n):
  rand_list = []
  from random import randint
  for a in range(a):
      for seq in DNA:
          r = randint(0, t-k)
          r_kmer = seq[r:r+k]
          rand_list.append(r_kmer)
  print(rand_list)

DNA = ['ACGAT', 'TAGCT', 'GAGCC']
t, k, n = 5, 3, 5 
random_kmers(DNA, t, k, n)

# Output -> ['GAT', 'AGC', 'GAG', 
             'ACG', 'GCT', 'GAG', 
             'CGA', 'AGC', 'AGC', 
             'GAT', 'GCT', 'AGC', 
             'CGA', 'AGC', 'GCC']
