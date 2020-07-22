        
def hamming_distance(pattern, seq):
    return sum(pattern_base1 != pattern_base2 for pattern_base1, pattern_base2 in zip(pattern, seq))

def total_minimum_distances(DNA, pattern, k):
    total_distance = 0
    for text in DNA:
        distances = []
        sub_string = []
        for i in range(len(text)):
            sub_string.append(text[i:i+k])
        for j in sub_string:
            distances.append(hamming_distance(pattern, j))
        total_distance = total_distance + min(distances)
    return total_distance        
    
pattern = 'TAA'
DNA = ['TTTATTT', 'CCTACAC', 'GGTAGAG' ]

total_minimum_distances(DNA, pattern, 3 )
# Output -> 3
