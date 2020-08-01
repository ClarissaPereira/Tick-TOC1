# Tick-TOC1
A repository of Python solutions to the code challenges from Unit 2 of Bioinformatics Algorithms (Stepik).

Before opening the code files, I recommend quickly reading about the [biology background](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/Biology%20Notes.md) of this project to get a clearer sense of the purpose of each program/function. I also highly recommend reading the [Algorithm Notes](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/Algorithm%20Notes.md) which will help you make a more informed decision on which of the functions in this repo are best suited to your purposes. 

The following Python files in this repo all serve similar purposes of identifying regulatory motifs in the upstream regions of a set of genes:
  * Brute Force Algorithms:
    * [basic motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/basic_motif_search.py)
    * [median motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/median%20motif%20finder.py)
  * Greedy Algorithms:
    * [greedy motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/greedy%20motif%20finder.py)
  * Randomised Algorithms:
    * [randomised motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/randomised%20motif%20search.py)
    * [Gibb's sampler](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/gibbs%20sampler.py)
    
   
The remaining files are smaller functions which contribute to the search algorithms.

The upsteam regions of ten dormancy-related genes from the Mycobacterium tuberculosis genome are included as test data for the randomised motif search and Gibb's sampler programs.
