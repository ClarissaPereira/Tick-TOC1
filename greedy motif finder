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

def generate_profile(motifs):
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
    import numpy
    kmer_probs = []
    for probs in range(len(profile)):
        a = numpy.prod(profile[probs])
        kmer_probs.append(a)
    return(kmer_probs)
