# A Monte Carlo algorithm which begins with a collection of randomly chosen k-mers Motifs (one from each sequence in DNA) 
# it then constructs a Position Probability Matrix for the motif collection and uses the PPM to generate new, more likely motifs
# This algorithm can be iterated N number of times to keep improving the motif collection 
# ultimately, it identifies the regulatory motif in each DNA sequence

import random 

def random_kmers(DNA, k):
# generates random k-mers found within DNA
    rand_list = []
    for sequence in DNA:
        l = len(sequence)
        r = random.randint(0, (l-k))
        rand_kmer = sequence[r:r+k]
        rand_list.append(rand_kmer)
    return (rand_list)

def most_common(motif_matrix):
# finds the most common nucleotide at any given position in a list of k-mers
    from collections import Counter
    return [Counter(col).most_common(1)[0][0] for col in zip(*motif_matrix)]

def consensus(motifs):
# uses most_common to create a k-mer motif of the most frequently occuring bases 
    consensus = []
    consensus.append(most_common(motifs))
    consensus_string = ''.join(consensus[0])
    return(consensus_string)

def hamming_distance(most_prob_kmer, kmer):
# calculates the number of mismatches between two k-mers
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(most_prob_kmer, kmer))

def calculate_score(motifs):
# calculates the total hamming distance between the consensus string and each k-mers in a list of motifs
# i.e. the lower the score, the closer the motif list is to the consensus string
    consensus_str = consensus(motifs)
    score = 0
    for motif in motifs:
        score = score + hamming_distance(consensus_str, motif)
    return score

def generate_profile(motifs, k):
# applies Laplace's rule of succession
# creates a position frequency matrix (PFM) based on a list of motif k-mers 
# creates a position probability matrix from the PFM
    import numpy
    profile = []
    probA, probC, probG, probT = [], [], [], []
           
    for base in range(k):
        countA, countC, countG, countT = 1, 1, 1, 1
        for motif in motifs:
            if motif[base] == "A":
                countA += 1
            elif motif[base] == "C":
                countC += 1
            elif motif[base] == "G":
                countG += 1
            elif motif[base] == "T":
                countT += 1
        totalcount = countA + countC + countG + countT
        probA.append(countA/totalcount)
        probC.append(countC/totalcount)
        probG.append(countG/totalcount)
        probT.append(countT/totalcount)
    profile.append(probA)
    profile.append(probC)
    profile.append(probG)
    profile.append(probT)
    t_profile = numpy.transpose(profile)
    return t_profile

def probability (profile, kmer):
# uses the position probability matrix of an existing k-mer motif list 
# calculates the probability of a possible kmer 
    product_probability = 1
    for base in range(len(kmer)):
        if kmer[base] == 'A':
            p = profile[base][0]
        elif kmer[base] == 'C':
            p = profile[base][1]
        elif kmer[base] == 'G':
            p = profile[base][2]
        elif kmer[base] == 'T':
            p = profile[base][3]
        product_probability = product_probability * p
    return(product_probability)   

def generate_motif(profile, DNA, k):
# uses the PPM to identify the k-mer within each DNA sequence with the highest product probability
    motif_list = []
    for sequence in DNA:
        prob_dict = {}
        for j in range(len(sequence)-k):
            kmer = sequence[j:j+k]
            prob_dict[kmer] = probability(profile, kmer)
        max_prob_kmer = max(prob_dict, key=prob_dict.get)
        motif_list.append(max_prob_kmer)
    return(motif_list)
            
    
def Randomised_Motif_Search(DNA, k):
# generates a profile based on existing motif collection and uses profile to generate a new motif collection
# if the new motif collection has a lower score that the existing one, it replaces it
# process continues until score fails to improve
    initial_motifs = random_kmers(DNA, k)
    BestMotifs = initial_motifs
    while True:
        Profile = generate_profile(BestMotifs, k)
        Motifs = generate_motif(Profile, DNA, k)
        if calculate_score(Motifs) < calculate_score(BestMotifs):
            BestMotifs = Motifs
        else:
            return BestMotifs
        
DNA=['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA'
     'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG'
     'TAGTACCGAGACCGAAAGAAGTATACAGGCGT'
     'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC'
     'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
k = 8
N = 1000
# number of iterations of the randomised algorithm
BestMotifs = Randomised_Motif_Search(DNA, k)
for n in range (1, N):
    Motifs = Randomised_Motif_Search(DNA, k)
    if calculate_score(Motifs) < calculate_score(BestMotifs):
            BestMotifs = Motifs
print(BestMotifs)
#print(calculate_score(BestMotifs))

# Output -> ['TCTCGGGG', 'CCAAGGTG', TACAGGCG', 'TTCAGGTG', 'TCCACGTG']
