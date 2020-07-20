import itertools
string = 'AGCT'
kmers = itertools.product(string, repeat = 3 )
new_list = (''.join(w) for w in list(kmers))
allk = list(new_list)

dna_strings = ['ATTTGGC','TGCCTTA','CGGTATC','GAAAATT']
patterns1 = []
for i in dna_strings:
    for j in range(len(i)-3+1):
        patterns1.append(i[j:j+3])

for i in patterns1:
    for j in allk:
        hamming = 0
        for l in range(len(j)):
            if i[l] != j[l]:
                hamming += 1
        if hamming <= 1:
            patterns1.append(j)

print(patterns1) 
