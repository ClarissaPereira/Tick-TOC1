# Regulatory Motif Search Algorithms
## Brute Force Algorithms 
These algorithms conduct an exhaustive search across a collection of DNA sequences to find all k-mer motifs with up to *d* mismatches. An example is motif enumeration which generates a neighborhood of k-mers with up to d-mismatches for each k-mer in the DNA sequences and then cycles through each k-mer neighborhood again to check for frequently repeated k-mers. Whilst guaranteed to return the correct motifs, brute force algorithms are very inefficient for large datasets and they also rely on the bioinformatician already knowing how much variation can be tolerated.

A major issue with the brute force algorithm is that biological data, whether in DNA arrays or computational arrays, is inherently noisy. In a noisy dataset, not all the DNA sequences will contain a regulatory motif which will cause the motif enumeration algorithm to fail as it relies on motif matches being present in each sequence. Direct comparison of k-mers using the Hamming Distance function also isn't possible due to the high levels of variation in motifs; for example, there might be only four mismatches between the regulatory motif sequence and two mutated motif k-mers found in separate sequences but eight mismatches between the two k-mers themselves. To compensate for this, the motif enumeration algorithm can be modified to find k-mers that match an ideal motif (instead of searching for k-mers that match each other). Since the ideal motif is unknown, the frequency of bases in a set of k-mers can be used to devise a **consensus string** that is most likely to be the un-mutated ideal motif. For position *i* in each k-mer, the most frequent base is selected and used to build up the consensus string. For instance, a list of k-mers *ACT*, *AGT*, *GGG* would produce the consensus string *AGT*.

This can be done by creating a **Position Frequency Matrix (PFM)** where each column represents position *i* of the collection of motif k-mers taken from different DNA sequences and the rows represent the frequency of each base in position *i*. For instance, this diagram shows how a position frequency matrix is created from the ten Drosophila NF-kB binding sites (the score is the total number of mismatches between each motif and the new consensus string):

<img src="https://user-images.githubusercontent.com/68158694/88489227-fb7e0d00-cf8a-11ea-88c0-6a159ec7f1f6.png" width=650 align=center>


PFMs also need to take into account that, in certain positions in the motif, multiple bases may have the same ability to bind to a transcription factor. In the example above, at position 7, both *C* and *T* demonstrate an equal ability to bind to the transcription factor. Therefore, a consensus string that reflects variation for the NF-kB binding site would be *TCGGGGA[C/T]TT[A/C/T]C*. To account for such variation in the consensus string, total entropy can be used in place of a simple score. Entropy is a measure of the uncertainty of a probability distribution. For each column in the PFM, the frequency of a particular base is converted into a probability (creating a **Position Probability Matrix or PPM**). Column entropy = the sum of (probability x log<sub>2</sub>probability) for each base in the column. All the column entropies are then added together to calculate the total entropy of the motif matrix. The lower the total entropy, the more conserved the motif matrix is and the better each motif matches the consensus string.

To improve efficiency of the brute force algorithm, we can instead search for a **median motif** which minimise the Hamming distance between itself and each DNA sequence k-mer. Instead of exploring all possible motifs in the collection of DNA sequences and then deriving the consensus string, this new algorithm will explore all potential k-mer consensus strings first and then find the best possible collection Motifs for each consensus string. Therefore, it will run fewer cycles through each DNA sequence, identifying the regulatory motif faster. However, for longer 20-mer motifs, this median motif finder algorithm is still not efficient enough. 

#### brute force algorithms:
* [motif enumeration](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/basic_motif_search.py)
* [median motif finder](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/median%20motif%20finder.py)

## Greedy Algorithms
Unlike a brute force algorithm which works iteratively through evey possible k-mer to find the regulatory motif, greedy algorithms work heuristically by choosing the optimal k-mer at each step as it attempts to find the overall optimal motif. For instance, a greedy algorithm would scan through the first DNA sequence, choose the most probable motif and add it to the motif matrix. This motif matrix is then used to scan through the second DNA sequence and choose the most probable motif based on the PPM of the motifs previously added to the matrix. Instead of starting out with a motif matrix and comparing all k-mers to it and repeating this process until a minimum Hamming Distance is found, the greedy algorithm builds up a motif matrix as it goes along. 

One drawback of using a greedy algorithm to search for motifs is that it trades accuracy for speed by not considering all possible options and never reconsidering a previous decision.

<img src="https://user-images.githubusercontent.com/68158694/88552381-ef8e5b80-d01b-11ea-973b-ea47aed6839d.png" width=250 align=left>

This limitation can be illustrated by applying a greedy algorithm to a simple graph search; the algorithm below traverses the graph and adds the number on each node to a running total. It aims to find a path that will produce the highest total. Just from observation, we know that the path with the highest total must include the *99* node because it's the largest number by far. So the ideal path would be (7 -> 3 -> 1 -> 99).

But by selecting the highest option at each step, never reconsidering a previous decision, and not evaluating all possible paths, the greedy algorith will instead choose the path (7 -> 12 -> 6 -> 9). The solution is found much faster but it is definitely not the overall optimal solution.

.    .    .    .   .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    <img src="https://user-images.githubusercontent.com/68158694/88552363-eac9a780-d01b-11ea-9027-39d064cf8b1a.gif" width=450 align=center>    .    .    .    .   .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .

To prevent an initial misstep derailing the entire algorithm, the first motif in the matrix must be chosen with care. First, the greedy algorithm does an exhaustive search on the first two strands of DNA to determine the best motif in these two strands. This motif becomes a seed node which influences the rest of the greedy search. Next, the algorithm searches each of the remaining DNA *N-2* strands for a motif which best matches the seed and the motifs that have already been found. 

Another method of improving algorithm accuracy involves the Cromwell Rule probability theory. In any observed data set, there is the possibility, especially with low-probability events or small data sets, that an event with nonzero probability does not occur. Its observed frequency is therefore zero; however, setting the empirical probability of the event equal to zero represents an inaccurate oversimplification that may cause problems. For instance, if our motifs matrix was [*AGT*, *ACT*], the consensus string would be *A[G/C]T* and the corresponding PPM would be:

[[1.0, 0.0, 0.0, 0.0], 

[0.0, 0.5, 0.5, 0.0],

[0.0, 0.0, 0.0, 1.0]]. 

Let's consider the probabilities of two possible motif k-mers:  *GTA* and *ATT* 
* P(*GTA*) = 0.0 x 0.0 x 0.0 = 0
* P(*ATT*) = 1.0 x 0.0 x 1.0 = 0

Both k-mers would have a probability of 0 despite ATT being only one nucleotide off from the consensus string and much closer to the ideal motif than *GTA*. By artificially adjusting the probability of rare events, these problems can be mitigated. We can do this by applying Laplace's Law of Succession and giving each nucleotide in the matrix a pseudocount which starts at 1 - altering the PFM. The new resulting PPM would look like:

[[1/2, 1/6, 1/6, 1/6], 

[1/6, 1/3, 1/3, 1/6],

[1/6, 1/6, 1/6, 1/2]]

And we can again compare the probabilities of two possible motif k-mers:  *GTA* and *ATT* 
* P(*GTA*) = 1/6 x 1/6 x 1/6 = 1/216
* P(*ATT*) = 1/2 x 1/6 x 1/2 = 1/24

This time, the algorithm does not discount *ATT* since, although both k-mers are not perfect matches, *ATT* has a much higher probability than *GTA*.

#### greedy algorithms:
* [greedy motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/greedy%20motif%20finder.py)
