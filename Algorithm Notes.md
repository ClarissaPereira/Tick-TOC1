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

The algorithm no longer discounts *ATT* since, although both k-mers are not perfect matches, *ATT* has a much higher probability than *GTA*.

#### greedy algorithms:
* [greedy motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/greedy%20motif%20finder.py)

## Randomised Algorithms
These algorithms employ a degree of randomness as part of their logic - for instance, they might generate random inputs and then perform deterministic computations on them. Randomised algorithms fall into one of two groups: Monte Carlo algorithms and Las Vegas algorithms. Monte Carlo algorithms have a fixed runtime (at the slight expense of accuracy) whereas Las Vegas algorithms always return a correct solution although runtime varies.

Randomness in a motif search seems counter-intuitive but it can be extremely useful for handling large datasets; whereas brute force algorithms must search through each possible option before returning a solution, randomised algorithms can reach a solution much faster. For instance, imagine we are looking for the number 7 in this list: [2, 14, 5, 9, 23, 42, 13, 7]. The order of this list is actually the worst-case scenario for this search because 7 is right at the end. 

A brute force algorithm for finding 7 would look like this:

```python
for i in range(len(number_list)):
  if number_list[i] == 7:
    return(7)
```

It would require eight iterations to find 7 whereas a randomised algorithm could find 7 more efficiently by instead randomly stumbling upon it:

```python
r = random.randint(0, len(number_list)-1)
while True:
  if number_list[r] == 7:
    return(7)
```

This randomised algorithm could require just one iteration to find 7. At worst, it would still require eight iterations- making it a more ideal option for searching than a brute force algorithm. 

A classic use of a Monte Carlo algorithm is to approximate the value of pi. The algorithm models a circle nested in a square of the same 'radius' so that the ratio of areas (circle/square) = pi/4. Points are randomly scattered within the square; if they fall inside the circle they turn green, otherwise they turn red. The number of green points divided by the total number of points is equivalent to the ratio of areas. Pi can then be calculated as (total # of green points / total # of points) x 4. Below, I have illustrated this concept using Python Turtle:

![ezgif com-video-to-gif](https://user-images.githubusercontent.com/68158694/89105283-bab54680-d417-11ea-9713-3f0a3a133644.gif)
  
This simple example highlights an important consideration for Monte Carlo algorithms; the algorithm must be iterated many times. Since the number of randomly generated points represents the total areas, estimation improves as more points are placed. Using 3,000 points, as this simulation has, gives a pi value of 3.137 (only accurate to 3sf). Similarly, for motif finding, a single iteration of a randomised search would return poor results; the algorithm must be repeated several hundreds of times to find the best approximate solution. 

In each iteration, *Randomised Motif Search* randomly selects one k-mer from each sequence in a list of DNA sequences and generates a PPM for the random motifs. This PPM is then used to identify one k-mer in each DNA sequence with the highest positional probability. These become the new motifs. If the new motifs are more conserved than the previous set of motifs, then it replaces them as the best set of motifs. The algorithm repeats until the new set of motifs are no longer an improvement on the best motifs. After 1,000 iterations of this process, the regulatory motifs in the DNA sequences are identified.

One drawback of *Randomised Motif Search* is that, in each iteration, all existing motifs are replaced which means that some correct motifs may be discarded and will have to be found again. A more cautious **Gibb's Sampler** algorithm, which replaces one motif at a time, could reduce runtime and increase the accuracy of the motif search. A Gibb's sampler randomly selects a single motif to be removed from a best motif set. It then uses a slightly more advanced random number generator to select a new motif from the DNA sequence (based on the PPM of the four remaining motifs). In essence, the randomness of the algorithm acts as a biased die - improving the selection of new motifs.

#### randomised algorithms:
  * [generate random k-mers](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/generate%20random%20k-mers.py)
  * [randomised motif search](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/randomised%20motif%20search.py)
  * [Gibb's sampler](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/gibbs%20sampler.py)

