# Tick-TOC1
A repository of Python solutions to the code challenges from Unit 2 of Bioinformatics Algorithms (Stepik).

Before opening the code files, I recommend quickly reading about the [biology background](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/Biology%20Notes.md) of this project to get a clearer sense of the purpose of each program/function. I also highly recommend reading the [Algorithm Notes](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/Algorithm%20Notes.md) which will help you make a more informed decision on which of the functions in this repo are best suited to your purposes. 

The upsteam regions of ten dormancy-related genes from the Mycobacterium tuberculosis genome are included as test data for the randomised motif search and Gibb's sampler programs. These two algorithms identified the following as potential regulatory motifs for the Dormancy Survival Regulator (DosR): 

'GACTTCAGGCCC', 'GACCCGCGGCCC', 'GACTTCCGGCGG', 'GAATCCCGGACC', 'GACCGACGTCCC', 'GACCTTCGGCCC', 'GACTTCTGTCCC', 'GACTTTCGGCCC', 'GACTAACGGCCC', 'GACCGCCTGGCC' 

These motifs have a low mismatch score of 24 indicating a high degree of conservation. The [consensus function](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/greedy%20motif%20finder.py) (taken from the greedy motif search file) identifies the median motif as **GACTTCCGGCCC**.

Comparison with [literature](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1992516/) reveals this 12-mer motif to be an accurate segment of the larger identified 20-mer motif TT[C/G]GG**GACT[T/A]TAGTCCC[G/C]** AA. The published, experimentally-identified motif and the motif identified by the algorithms in this repo have just four mismatches between them - a 75% success rate which could be improved with more data. 

The following Python files in this repo all serve similar purposes of identifying regulatory motifs in the upstream regions of a set of genes:
### Brute Force Algorithms:
  * [basic motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/basic_motif_search.py)
  * [median motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/median%20motif%20finder.py)
### Greedy Algorithms:
  * [greedy motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/greedy%20motif%20finder.py)
### Randomised Algorithms:
  * [randomised motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/randomised%20motif%20search.py)
  * [Gibb's sampler](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/gibbs%20sampler.py)
    
   
The remaining files are smaller functions which contribute to the search algorithms.


