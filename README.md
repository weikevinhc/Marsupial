# MarSuPial
**Estimating recombinant fraction, genetic distance, and recombination rate from Marker Selected Pools**

Marker selected pools are pooled-sequencing of marker selected offsprings from a two generation recombinant backcross scheme. The MarSuPial package takes filtered and polarized allele-specific count data and infers the recombinant fraction, genetic distance and recombination rates from the allele-frequency decay across the chromosome. MarSuPial does not genotype, filter, or polarize counts. We recommend following GATK best practices for genotyping.
MarSuPial also simulates read count data for allele-frequency decay provided the chromosome-wide recombination rate.  

## Installation

Download the package MarSuPial_1.0.tar.gz and install using command line:
    <p>$ R CMD INSTALL MarSuPial_1.0.tar.gz</p>

Then, in R, load the package:
    <p>> library(MarSuPial)</p>

OR Install in R by first extracting the tar ball and:
    <p>> install.packages(MarSuPial)</p>
    <p>> library(MarSuPial)</p>

## Usage
The MarSuPial package contains a collection of functions to analyzing, estimating, and simulating recombination in a two generation backcross. Refer to the MarSuPial_scripts folder for the code.

### Analytical relationships:
r2d.locus(rate, l, size, pos, start, end)<br>
Returns a list of the chromosome-wide windows and the recombination rate and genetic distance in cM from a locus (l). Rate (in cM/Mb) must be either a r function describing the relationship between recombination rate or chromosome position or a bed file (loaded as a dataframe) of equally sized chromosome windows. If pos (in Mb) is a single value, it is used as the window size of evenly spaced non-overlapping windows. If pos (in Mb) is a vector of numbers, they are used as positions to determine rate and distance. If a recombination rate map function is provided, the start and end positions of the recombination rate can also be inputted; positions < start and > end will have rates of 0 to emulate heterochromatic/centromeric/telomeric suppression.

d2D(d, method)<br>
Converts genetic distance (d) into recombinant fraction (D), using either the "haldane" or "kosambi" mapping functions (method).

D2d(D, method)<br>
Converts recombinant fraction (D) into genetic distance (d), using either the "haldane" or "kosambi" mapping functions (method).

D2AF(D, fitness, ploidy, mendel_rate)<br>
Converts recombinant fraction (D) to allele frequency (AF) after taking taking into account the fitness differential (fitness) of the selected locus, the ploidy, and the expected mendelian ratio (mendel_rate).

AF2D(AF, fitness, ploidy, mendel_rate)<br>
Converts allele frequency (AF) to recombinant fraction (D) after taking taking into account the fitness differential (fitness) of the selected locus, the ploidy, and the expected mendelian ratio (mendel_rate).

r2d.3pt(ratefunc, l1, l2, size, winsize = 0.01, start, end)<br>
Same as r2d.locus, but for two loci under selection. Rate must be a function.

d2AF.3pt(list3pt, fitness1, fitness2, method)<br>
Takes an r2d.3pt outlut (list), and generates the alleles frequency given the fitness differential of the two loci.

### Allele frequency inference:
AFfitsplines(sites, freq, depth, locus, locus_freq, df)<br>
Fits a smoothing curve over chromosome-wide allele frequency (freq) from read count data at informative positions. Positions (sites) are SNPs where the lines are homozygous for different nucleotides. The coverage/readcount (depth), and position of selection (locus), and expected selection frequency (locus_freq) must provided. The extent of smoothing is determine by the degree of freedom (df). By default, a cross validation method is used, but this may overfit the data. For a chromsome window size quivalent use, df = chromosome size/ window size - 4. Returns a list object (AFfit) which contains two splines fits, one on either side of the selected locus.

AFfitloess(sites, freq, depth, locus, locus_freq, span)<br>
Same as AFfitsplines, but uses loess smoothing. The smoothing factor is controlled by span. 

AFwinlm(sites, freq, window, slide)<br>
Estimates the allele frequency at sliding windows using linear regressions in each window. The size of the window is determeined by the window parameter and how far each window is apart by the slide parameter. Must be in Mb. 

predictAF(AFfit, sites)<br>
Takes an AFfit object (from AFfitsplines or AFfitloess), and generates the predicted allele frequency at any given position (sites).

predictm(AFfit, sites)<br>
Takes an AFfit object from AFfitsplines, and generates the slope of the allele frequency at any given position (sites).

d2rslope(d, sites, locus)<br>
Converts genetic distance to recombination rate by taking the slope of the former in discrete intervals. Regions to the left and right of the selected locus is expected to have negative and positive slopes. So slopes of regions left of the locus are multiplied by -1. 

### Crossover simulations of two generation recombinant backcross:
COsim(r2d, fintess, n, trials)<br>
Takes the output of r2d.locus, and simulates crossing over for n number of chromosomes and repeated trials number of times, given different fitness differential (fitness must be between -1 and 1). Produces a matrix of allele frequencies where each column is a site in the r2d list, and each rwo is a trial. Allele frequency here refers to true allele frequency in the pool. This simulation method assumes no crossover interference.

COsim.ind(r2d, fitness, n)<br>
Instead of generating the chromosome-wide allele frequency of a pool, this function generates a matrix of the genotypes of individual chromosomes where each row is an individual and each row a position (True and False for the two allele states). This is equivalent to using COsim with n = 1, but is faster. However it takes a long time to generate large number of trials.

AF2WGSsim(AF, dp)<br>
Takes the vector of allele frequency and simulates read counts for each. Returns a list object with counts of the allele and depth of the site. Depth is randomly sampled using a Poisson process with mean of dp. The allele frequency is sampled with binomial. 

D2WGSsim(D, fitness, n, dp, mendel_rate)<br>
A rough and fast way to simulate allele frequency from read count data given recombinant fraction (D) and fitnesss differential. For each D, a binomial sampling process is used to determine the number of alleles in the pool, followed by sequnece depth smapling by poisson, and then another binomial for the read count of the allele. This should only be applied when n is large (>1000) as it samples the allele frequency at each site independently. 

### Examples:
See the vignette (Rmark down or PDF).
