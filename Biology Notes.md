# Tick-TOC1 - Identifying Molecular Clocks
## Background
The behaviour and activity of animals, plants, and even bacteria is controlled by an internal timekeeper known as a **circadian clock**. This biochemical oscilllator has a steady rhythm and synchronises with fluctuating light levels. Even in the complete absence of light, our circadian clocks maintain an approximate 24 hour cycle. 

Identifying molecular clocks can help us better understand inherited disorders, such as **delayed sleep-phase syndrome**. The activity of these miniscule timepieces could also explain why heart attacks occur more often in the morning whereas asthma attacks are more common at night.

Given that so much of their activity is centred around solar time, molecular clocks are even more vital to the survival of plants. Each plant cell tracks time independently with just three regulatory genes (**LHY**, **CCA1**, and **TOC1**) involved in a negative-feedback loop. During the day, sunlight activates LHY and CCA1 transcription; in turn, these genes repress transcription of TOC1. As night falls, light levels dip and LHY and CCA1 are no longer upregulated. In the absence of repression, TOC1 is transcribed and promotes expression of LHY and CCA1. Each of these regulatory genes encode **transcription factors** that control expresion of many other genes by binding to **regulatory motifs** in their upstream regions.  

## Identification Techniques
A wet lab approach to identifying regulatory motifs is to create a DNA array (a collection of short DNA probe molecules attached to a silicone chip). Fluorecently labelled slices of the investigated genome are added to the array and bind to probes via complementary base pairing. The greater the expression of a gene, the greater the intensity of the fluorescent signal. Since an array may contain millions of probes, biologists can measure the expression of many genes in a single array experiment. For instance, in 2000, chronobiologist Steve Kay identified clock genes in thale cress (*Arabidopsis thaliana*) with a DNA array experiment that measured the expression of 8,000 genes. Other regulatory motifs, such as the **NF-kB binding site** in Drosophila genomes, are not as well-conserved and therefore are far more elusive in a DNA sequence. Computational arrays can be used to identify such motifs. 

One difficulty with using a computational method is that, unlike DnaA boxes which cluster together into clumps, regulatory motifs are scattered throughout the genome. Motifs only occur at least once (with variations and mutations) in many different genoke regions.  As a result, the [algorithm designed to locate DnaA boxes](https://github.com/ClarissaPereira/Finding-Ori/blob/master/Final%20DnaA%20Box%20Finder.py) in the origin of replication is inadequate at modelling the biological problem of regulatory motif finding. 

## Brute Force Algorithms 
These algorithms conduct an exhaustive search across a collection of DNA sequences to find all k-mer motifs with up to *d* mismatches. An example is motif enumeration which generates a neighborhood of k-mers with up to d-mismatches for each k-mer in the DNA sequences and then cycles through each k-mer neighborhood again to check for frequently repeated k-mers. Whilst guaranteed to return the correct motifs, brute force algorithms are very inefficient for large datasets and they also rely on the bioinformatician already knowing how much variation can be tolerated.

Direct comparison of k-mers using the Hamming Distance function also isn't possible due to the high levels of variation in motifs; for example, there might be only four mismatches between the regulatory motif sequence and two mutated motif k-mers found in separate sequences but eight mismatches between the two k-mers themselves. To compensate for this, the motif enumeration algorithm can be modified to find k-mers that match an ideal motif (instead of searching for k-mers that match each other). Since the ideal motif is unknown, the frequency of bases in a set of k-mers can be used to devise a **consensus string** that is most likely to be the un-mutated ideal motif. For position *i* in each k-mer, the most frequent base is selected and used to build up the consensus string. For instance, a list of k-mers *ACT*, *AGT*, *GGG* would produce the consensus string *AGT*.

This can be done by creating a **Position Frequency Matrix (PFM)** where each column represents position *i* of the collection of motif k-mers taken from different DNA sequences and the rows represent the frequency of each base in position *i*. For instance, this diagram shows how a position frequency matrix is created from the ten Drosophila NF-kB binding sites (the score is the total number of mismatches between each motif and the new consensus string):

![motif snip](https://user-images.githubusercontent.com/68158694/88489227-fb7e0d00-cf8a-11ea-88c0-6a159ec7f1f6.png)

PFMs also need to take into account that, in certain positions in the motif, multiple bases may have the same ability to bind to a transcription factor. In the example above, at position 7, both *C* and *T* demonstrate an equal ability to bind to the transcription factor. Therefore, a consensus string that reflects variation for the NF-kB binding site would be *TCGGGGA[C/T]TT[A/C/T]C*. To account for such variation in the consensus string, total entropy can be used in place of a simple score. 

Entropy is a measure of the uncertainty of a probability distribution. For each column in the PFM, the frequency of a particular base is converted into a probability (creating a **Position Probability Matrix**). Column entropy = the sum of (probability x log<sub>2</sub>probability) for each base in the column. All the column entropies are then added together to calculate the total entropy of the motif matrix. 

A major issue with such an algorithm is that biological data, whether in DNA arrays or computational arrays, is inherently noisy. In a noisy dataset, not all the DNA sequences will contain a regulatory motif which will cause the motif enumeration algorithm to fail as it relies on motif matches being present in each sequence. 

### brute force algorithms:
* [motif enumeration (basic motif search)](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/basic_motif_search.py)



