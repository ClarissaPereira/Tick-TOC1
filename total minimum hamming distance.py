def sub_strings(text, k):
    for i in range(len(text) - k + 1):
        yield text[i:i+k]
        
def hamming_distance(pattern, seq):
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(pattern, seq))

def total_minimum_distances(DNA, pattern, k):
    distance = 0
    for text in DNA:
        for sub_string in sub_strings(text, k):
            min_hamming = float('inf')
            if min_hamming > hamming_distance(pattern, sub_string):
                min_hamming = hamming_distance(pattern, sub_string)
        distance = distance + min_hamming
    return distance        
    
pattern = 'TAA'
DNA = ['TTTATTT', 'CCTACAC', 'GGTAGAG' ]

total_minimum_distances(DNA, pattern, 3 )
