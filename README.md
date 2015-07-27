# Salix

Gene-dropping of sequence variants within pedigrees

Brief description of the files:
* *cm_length.txt* contains the length of each chromosome in centiMorgans.  A Morgan is a measure of genetic distance which represents the expected number of crossovers between two markers.
* *run.sh* is an example Bash script for the entire process.

## Map distance and recombination fractions

If markers are spaced at 0.1 CentiMorgans (cM) this signifies that there's a ```0.1/100=.001``` probability of adjacent markers recombining according to the Morgan map function. Note that the probabilities are considered independent.

The following is taken verbatim from [Wallace 2003](https://www-gene.cimr.cam.ac.uk/staff/wallace/pub/thesis.pdf).

### Morgan
The simplest map function, known as the Morgan map function, assumes chromosomal segments can have at most one crossover, and that the probability of a crossover is proportional to the length of the segment. The probability of a crossover in m map units is therefore 2m and
```
θ= (1−p0)/2 = (1−(1−2m)) = m
```
Note that the function is only valid for m < 1/2 (otherwise we would have θ > 1/2) so is not applicable for long segments of chromosome.

### Haldane
The Haldane map function assumes crossovers occur at random, independently of one another. This implies a Poisson process, with rate 2m in a segment of length m.
So ```p0 = exp(−2m)``` and
```
θ= (1−exp(−2m)) / 2
```
However, observations show that crossovers do not tend to occur completely independently.

### Kosambi
Due to interference, the probability of having two crossovers very close to each other is less than that predicted by the Haldane function, and more complex arguments that take account of this lead to the Kosambi function
```
θ= (exp(4m)−1) / 2(exp(4m) + 1)
```
Other functions have also been proposed, but the above three are the sim- plest. Figure 2.2 shows that all functions produce very similar results for the small distances (< 10cM) which are generally used when conducting linkage analyses, as described in section 2.4.

| Chr | Centimorgans | Morgans |
| --- | ---| --- |
| 1 | 286.279234 |  2.86279234 |
| 2 | 268.839622 |  2.68839622 |
| 3 | 223.361095 |  2.23361095 |
| 4 | 214.688476 |  2.14688476 |
| 5 | 204.089357 |  2.04089357 |
| 6 | 192.039918 |  1.92039918 |
| 7 | 187.220500 |  1.872205 |
| 8 | 168.003442 |  1.68003442 |
| 9 | 166.359329 |  1.66359329 |
| 10 | 181.144008 | 1.81144008 |
| 11 | 158.218650 | 1.5821865 |
| 12 | 174.679023 | 1.74679023 |
| 13 | 125.706316 | 1.25706316 |
| 14 | 120.202583 | 1.20202583 |
| 15 | 141.860238 | 1.41860238 |
| 16 | 134.037726 | 1.34037726 |
| 17 | 128.490529 | 1.28490529 |
| 18 | 117.708923 | 1.17708923 |
| 19 | 107.733846 | 1.07733846 |
| 20 | 108.266934 | 1.08266934 |
| 21 | 62.786478 |  0.62786478 |
| 22 | 74.109562 |  0.74109562 |


PhD: "Genetic Susceptibility to leprosy: methodological issues in a linkage analysis of extended pedigrees from Karonga District, Malawi", LSHTM, 2003 [thesis](https://www-gene.cimr.cam.ac.uk/staff/wallace/pub/thesis.pdf)

[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/APLevine/Salix?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
