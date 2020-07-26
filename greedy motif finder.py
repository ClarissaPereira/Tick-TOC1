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

def generate_profile(motifs):
  # creates a position frequency matrix (PFM) based on a list of motif k-mers 
  # creates a position probability matrix from the PFM
    import numpy
    profile = []
    probA, probC, probG, probT = [], [], [], []
    for j in range(len(motifs[0])):
        countA, countC, countG, countT = 0, 0, 0, 0
        for motif in motifs:
            if motif[j] == "A":
                countA += 1
            elif motif[j] == "C":
                countC += 1
            elif motif[j] == "G":
                countG += 1
            elif motif[j] == "T":
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

def probability (profile):
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

def greedy_motif_search(DNA, k):
  # a greedy algorithm which returns a list of best motifs (one from each DNA sequence)
  # best motifs are k-mers which most closely match each other
    
    # initialise a few variables
    for sequence in DNA:
        best_motif.append(sequence[:k])
    first_seq_kmers = []
    first_seq = DNA[0]
    for j in range(len(first_seq) - k + 1):
        first_seq_kmers.append(first_seq[j:j+k])
        
    # greedy algorithm begins here
    for x in first_seq_kmers:
        kmer_list = []
        kmer_list.append(x)
        for string in DNA[1:len(DNA)]:
            prob_dict = {}
            for n in range(len(string) - k + 1):
                kmer = (string[n:n+k])
                motifs_profile = generate_profile(kmer_list)
                prob_dict[kmer] = probability(motifs_profile, kmer)
            #print(prob_dict)
            if all(value == 0 for value in prob_dict.values()):
                kmer_list.append(string[:k])
            else:
                max_prob_kmer = max_key = max(prob_dict, key=prob_dict.get)
                kmer_list.append(max_prob_kmer)
        print(kmer_list)
        print(calculate_score(kmer_list))
        if calculate_score(kmer_list) < calculate_score(best_motif):
            best_motif = kmer_list
    print(best_motif)

greedy_motif_search(DNA = ['GGCGTTCAGGCA', 
                           'AAGAATCAGTCA', 
                           'CAAGGAGTTCGC', 
                           'CACGTCAATCAC', 
                           'CAATAATATTCG'], k = 3)

# Output -> ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']





