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

### brute force algorithms:
* [motif enumeration (basic motif search)](https://github.com/ClarissaPereira/Tick-TOC1/blob/master/basic_motif_search.py)



