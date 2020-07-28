import random 

def random_kmers(DNA, k):
    rand_list = []
    for sequence in DNA:
        l = len(sequence)
        r = random.randint(0, (l-k))
        rand_kmer = sequence[r:r+k]
        rand_list.append(rand_kmer)
    return (rand_list)

def most_common(motif_matrix):    
    from collections import Counter
    return [Counter(col).most_common(1)[0][0] for col in zip(*motif_matrix)]

def consensus(motifs):
    consensus = []
    consensus.append(most_common(motifs))
    consensus_string = ''.join(consensus[0])
    return(consensus_string)

def hamming_distance(most_prob_kmer, kmer):
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(most_prob_kmer, kmer))

def calculate_score(motifs):
    consensus_str = consensus(motifs)
    score = 0
    for motif in motifs:
        score = score + hamming_distance(consensus_str, motif)
    return score

def generate_profile(motifs, k):
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
    initial_motifs = random_kmers(DNA, k)
    BestMotifs = initial_motifs
    while True:
        Profile = generate_profile(BestMotifs, k)
        Motifs = generate_motif(Profile, DNA, k)
        if calculate_score(Motifs) < calculate_score(BestMotifs):
            BestMotifs = Motifs
        else:
            return BestMotifs
        
DNA=[]
k = 
N = 1000
BestMotifs = Randomised_Motif_Search(DNA, k)
for n in range (1, N):
    Motifs = Randomised_Motif_Search(DNA, k)
    if calculate_score(Motifs) < calculate_score(BestMotifs):
            BestMotifs = Motifs
print(BestMotifs)
print(calculate_score(BestMotifs))
